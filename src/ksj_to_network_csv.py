#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KSJ道路中心線GeoJSON → 道路ネットワークデータ変換スクリプト

対象範囲・速度設定を柔軟に指定できる統合版。
単一1次メッシュ / 複数1次メッシュ / 全国 を同一スクリプトで処理する。

【使い方】
  # 単一メッシュ（市区町村・都道府県詳細）
  python3 ksj_to_network_csv.py --mesh 5238 --case 5238_shizuoka
  python3 ksj_to_network_csv.py --mesh 5238 --case 5238_shizuoka_walk --mode walk

  # 複数メッシュ（地方ブロック: フィルターあり推奨）
  python3 ksj_to_network_csv.py \\
    --meshes 5338,5339,5340,5439,5440 --case kanto_block --filter

  # 複数メッシュ（都道府県: フィルターなし）
  python3 ksj_to_network_csv.py \\
    --meshes 5137,5138,5237,5238,5337,5338 --case shizuoka_pref

  # 全国（フィルターあり推奨）
  python3 ksj_to_network_csv.py --nationwide --filter
  python3 ksj_to_network_csv.py --nationwide --case nationwide_nofilt

【対象範囲オプション（いずれか1つ必須）】
  --mesh CODE        単一1次メッシュコード（例: 5238）
  --meshes CODE,...  複数1次メッシュ（カンマ区切り、例: 5338,5339,5340）
  --nationwide       全国全メッシュ（geojson/ 以下を全件 glob）

【速度モード（--mode）】
  vehicle（デフォルト）: 道路種別別速度（表9-2 DID内基準）
    国道=35, 都道府県道=30, 市区町村道=20, 高速=80 km/h
  walk: 徒歩速度（--walk-kmh で指定、デフォルト3.6 km/h）一律

【フィルター（--filter）】
  N13_003 in [1,2,4]（国道・都道府県道・高速）or N13_006 in [3,4,5]（幅員5.5m以上）
  地方ブロック・全国ではメモリ削減のため --filter 推奨。都道府県以下はなし推奨。

【出力（01_MakeNetwork/{case}/ に格納）】
  KSJ_N13-24_{case}_道路リンク.csv      リンクCSV（双方向展開、DRM3003形式準拠）
  KSJ_N13-24_{case}_道路リンク.parquet  リンクGeoParquet（正方向のみ、EPSG:6668）
  KSJ_N13-24_{case}_道路ノード.csv      ノードCSV（node_id, lon, lat）
  KSJ_N13-24_{case}_道路ノード.parquet  ノードGeoParquet（EPSG:6668）

【処理フロー】
  [Pass 1] GeoJSONスキャン → （フィルター）→ ノードID付与
  [Pass 2] 再スキャン → リンクCSV出力 + Parquet用データ収集
  [後処理] リンクParquet / ノードCSV + ノードParquet 出力
"""

import argparse
import json
import math
import os
import time
from pathlib import Path

import geopandas as gpd
from shapely.geometry import LineString, Point

# ============================
# 引数解析
# ============================
parser = argparse.ArgumentParser(
    description='KSJ道路ネットワーク変換（統合版: 単一/複数/全国メッシュ対応）')

# ── 対象範囲（いずれか1つ必須）──
range_group = parser.add_mutually_exclusive_group(required=True)
range_group.add_argument('--mesh', metavar='CODE',
                         help='単一1次メッシュコード（例: 5238）')
range_group.add_argument('--meshes', metavar='CODE,...',
                         help='複数1次メッシュ（カンマ区切り、例: 5338,5339,5340）')
range_group.add_argument('--nationwide', action='store_true',
                         help='全国全メッシュ（geojson/ 以下を全件 glob）')

# ── その他オプション ──
parser.add_argument('--case', default=None,
                    help='ケース名（出力サブディレクトリ名）。'
                         '省略時: --mesh → メッシュコード, --nationwide → nationwide')
parser.add_argument('--filter', action='store_true',
                    help='フィルター適用: 都道府県道以上(N13_003 in 1,2,4) or 幅員5.5m以上(N13_006 in 3,4,5)')
parser.add_argument('--mode', default='vehicle', choices=['vehicle', 'walk'],
                    help='速度モード: vehicle=道路種別別速度（デフォルト）, walk=一律 --walk-kmh km/h')
parser.add_argument('--walk-kmh', type=float, default=3.6,
                    help='walk モード時の速度（km/h）。デフォルト: 3.6')
parser.add_argument('--city', default=None,
                    help='市区町村名でリンク・ノードParquetをクリップ（N03-20250101.geojsonを使用）'
                         '（--prefと組み合わせで都道府県絞り込み可）例: 高松市 / 焼津市')
parser.add_argument('--pref', default=None,
                    help='都道府県名でリンク・ノードParquetをクリップ（prefecture.parquetを使用）'
                         '例: 香川県 / 鳥取県')
parser.add_argument('--out-dir', default=None,
                    help='出力先ディレクトリ（省略時: 01_MakeNetwork/{case}/）')
args = parser.parse_args()

# ── 対象範囲の決定 ──
NATIONWIDE = args.nationwide
if args.meshes:
    MESH_LIST = [m.strip() for m in args.meshes.split(',')]
elif args.mesh:
    MESH_LIST = [args.mesh]
else:
    MESH_LIST = None   # nationwide: ファイル一覧から動的取得

# ── ケース名デフォルト ──
if args.case:
    CASE_NAME = args.case
elif NATIONWIDE:
    CASE_NAME = 'nationwide'
elif args.mesh:
    CASE_NAME = args.mesh
else:
    parser.error('--meshes 使用時は --case が必須です')

USE_FILTER = args.filter
MODE       = args.mode
WALK_KMH   = args.walk_kmh
CLIP_CITY  = args.city
CLIP_PREF  = args.pref

# ============================
# 設定
# ============================
BASE_DIR    = Path(__file__).parent.parent
GEOJSON_DIR = BASE_DIR / 'input'
OUT_DIR     = Path(args.out_dir) if args.out_dir else BASE_DIR / 'network' / CASE_NAME
OUT_DIR.mkdir(parents=True, exist_ok=True)

CRS = 'EPSG:6668'   # JGD2011 地理座標系

# クリップ用ポリゴンソース
CITY_PARQUET = BASE_DIR / 'data' / 'city.parquet'
PREF_PARQUET = BASE_DIR / 'data' / 'prefecture.parquet'

OUT_LINKS_CSV     = OUT_DIR / f'KSJ_N13-24_{CASE_NAME}_道路リンク.csv'
OUT_LINKS_PARQUET = OUT_DIR / f'KSJ_N13-24_{CASE_NAME}_道路リンク.parquet'
OUT_NODES_CSV     = OUT_DIR / f'KSJ_N13-24_{CASE_NAME}_道路ノード.csv'
OUT_NODES_PARQUET = OUT_DIR / f'KSJ_N13-24_{CASE_NAME}_道路ノード.parquet'

# フィルター条件（--filter 指定時のみ使用）
ROAD_CLASS_OK = {'1', '2', '4'}   # 国道・都道府県道・高速自動車国道等
WIDTH_OK      = {'3', '4', '5'}   # 幅員5.5m以上

# N13_003 → 速度（km/h）vehicle モード
SPEED_KMH = {'1': 35, '2': 30, '3': 20, '4': 80, '5': 20, '6': 20}

# N13_003 → DRM道路種別コード
ROAD_TYPE = {'1': 3, '2': 5, '3': 7, '4': 1, '5': 7, '6': 7}

# N13_003 → 自専道フラグ（1=自専道, 2=一般道）
JISEN = {'4': 1}

PROG_STEP = 500_000   # 進捗表示間隔（フィーチャ数）


# ============================
# ユーティリティ
# ============================
def haversine_m(lon1, lat1, lon2, lat2):
    R = 6_371_000
    d_lat = math.radians(lat2 - lat1)
    d_lon = math.radians(lon2 - lon1)
    a = (math.sin(d_lat / 2) ** 2
         + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2))
         * math.sin(d_lon / 2) ** 2)
    return 2 * R * math.asin(math.sqrt(a))


def line_length_m(coords):
    return sum(
        haversine_m(coords[i][0], coords[i][1],
                    coords[i+1][0], coords[i+1][1])
        for i in range(len(coords) - 1)
    )


def mesh2_code(lon, lat):
    """経緯度 → 2次メッシュコード（6桁）"""
    p = int(lat * 1.5)
    u = int(lon - 100)
    q = int((lat - p / 1.5) * 12)
    r = int((lon - (100 + u)) * 8)
    return f'{p:02d}{u:02d}{q}{r}'


def passes_filter(props):
    """フィルター条件を満たすか判定"""
    n13_003 = str(props.get('N13_003', ''))
    n13_006 = str(props.get('N13_006', ''))
    return n13_003 in ROAD_CLASS_OK or n13_006 in WIDTH_OK


def find_geojson_files(mesh_list):
    """1次メッシュコードリストに対応するGeoJSONファイルパスを返す"""
    files = []
    for mesh in mesh_list:
        found = list(GEOJSON_DIR.glob(f'*/N13-24_{mesh}.geojson'))
        if found:
            files.extend(sorted(found))
        else:
            print(f'  警告: メッシュ {mesh} のGeoJSONが見つかりません '
                  f'（検索パス: {GEOJSON_DIR}/*/N13-24_{mesh}.geojson）')
    return files


# ============================
# QML出力
# ============================
def write_road_links_qml(qml_path):
    """道路リンクParquet用 QGIS カテゴリシンボルQML（N13_003別ライン色）を生成する。"""
    _DP = ('<data_defined_properties><Option type="Map">'
           '<Option name="name" type="QString" value=""/>'
           '<Option name="properties"/>'
           '<Option name="type" type="QString" value="collection"/>'
           '</Option></data_defined_properties>')
    # (value, label, color, line_width, pass)
    # pass: 描画順制御（数値が大きいほど前面）高速=3, 国道=2, 都道府県道=1, それ以下=0
    categories = [
        ('4', '高速自動車国道等', '0,48,135,255',  '0.7', '3'),
        ('1', '国道',           '230,0,38,255',   '0.5', '2'),
        ('2', '都道府県道',      '26,115,232,255', '0.3', '1'),
        ('3', '市区町村道等',    '40,167,69,255',  '0.1', '0'),
        ('5', 'その他',          '136,136,136,255','0.1', '0'),
        ('6', '不明',            '170,170,170,255','0.1', '0'),
    ]
    cats_xml  = '\n'.join(
        f'      <category symbol="{i}" value="{v}" label="{lbl}" render="true"/>'
        for i, (v, lbl, _, _, _) in enumerate(categories)
    )
    syms_xml  = '\n'.join(
        f'      <symbol name="{i}" type="line" alpha="1" clip_to_extent="1" is_animated="0" frame_rate="10">\n'
        f'        {_DP}\n'
        f'        <layer class="SimpleLine" enabled="1" pass="{pas}" locked="0">\n'
        f'          <Option type="Map">\n'
        f'            <Option name="line_color" type="QString" value="{col}"/>\n'
        f'            <Option name="line_style"  type="QString" value="solid"/>\n'
        f'            <Option name="line_width"  type="QString" value="{wid}"/>\n'
        f'            <Option name="line_width_unit" type="QString" value="MM"/>\n'
        f'          </Option>\n'
        f'        </layer>\n'
        f'      </symbol>'
        for i, (_, _, col, wid, pas) in enumerate(categories)
    )
    qml = f"""<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" styleCategories="Symbology">
  <renderer-v2 type="categorizedSymbol" attr="N13_003" forceraster="0" symbollevels="1" usingSymbolLevels="1" enableorderby="0">
    <categories>
{cats_xml}
    </categories>
    <symbols>
{syms_xml}
    </symbols>
    <rotation/>
    <sizescale/>
  </renderer-v2>
  <blendMode>0</blendMode>
  <featureBlendMode>0</featureBlendMode>
  <layerOpacity>1</layerOpacity>
</qgis>
"""
    with open(qml_path, 'w', encoding='utf-8') as f:
        f.write(qml)


# ============================
# メイン処理
# ============================
def main():
    t0 = time.time()

    # ── GeoJSONファイル一覧 ────────────────────────────────────
    if NATIONWIDE:
        geojson_files = sorted(GEOJSON_DIR.glob('*/N13-24_*.geojson'))
        mesh_desc = f'全国（{len(geojson_files)}ファイル）'
    else:
        geojson_files = find_geojson_files(MESH_LIST)
        mesh_desc = f'{MESH_LIST}（{len(geojson_files)}ファイル）'

    total_files = len(geojson_files)

    print(f'対象範囲          : {mesh_desc}')
    print(f'出力先            : {OUT_DIR}')
    print(f'フィルター        : {"あり（都道府県道以上 or 幅員5.5m以上）" if USE_FILTER else "なし（全道路）"}')
    print(f'速度モード        : {MODE}' + (f'（{WALK_KMH} km/h）' if MODE == 'walk' else ''))

    if total_files == 0:
        print('エラー: 対象ファイルが見つかりません。')
        return

    # ── [1/4] Pass 1: ノードID付与 ─────────────────────────────
    print('\n[1/4] Pass 1: ノードID付与（全ファイルスキャン）...')
    coord_to_nid: dict[tuple, str] = {}
    mesh_counters: dict[str, int]  = {}
    pass1_total    = 0
    pass1_filtered = 0

    for fi, fp in enumerate(geojson_files, 1):
        with fp.open(encoding='utf-8') as f:
            data = json.load(f)
        features = data['features']
        pass1_total += len(features)

        for feat in features:
            if USE_FILTER and not passes_filter(feat['properties']):
                continue
            pass1_filtered += 1
            coords = feat['geometry']['coordinates']
            for pt in (coords[0], coords[-1]):
                key = (round(pt[0], 9), round(pt[1], 9))
                if key not in coord_to_nid:
                    mesh = mesh2_code(pt[0], pt[1])
                    mesh_counters[mesh] = mesh_counters.get(mesh, 0) + 1
                    coord_to_nid[key] = f'{mesh}{mesh_counters[mesh]:05d}'

        if fi % max(1, total_files // 20) == 0 or fi == total_files:
            print(f'  {fi:>3}/{total_files} files  '
                  f'総計: {pass1_total:>10,}  フィルター後: {pass1_filtered:>8,}  '
                  f'ノード数: {len(coord_to_nid):>8,}  ({time.time()-t0:.0f}s)')

    print(f'\n  フィルター後フィーチャ数 : {pass1_filtered:,}')
    print(f'  ユニークノード数         : {len(coord_to_nid):,}')
    print(f'  使用2次メッシュ数        : {len(mesh_counters):,}')
    print(f'  経過時間: {time.time()-t0:.0f}s')

    # ── [2/4] Pass 2: リンクCSV出力 + Parquet用データ収集 ──────
    print(f'\n[2/4] Pass 2: リンクCSV出力・Parquet用データ収集 ...')
    print(f'      CSV: {OUT_LINKS_CSV}')

    CSV_HEADER = (
        '2次ﾒｯｼｭ,ﾉｰﾄﾞ1,ﾉｰﾄﾞ2,管理者,道路種別,路線番号,市町村コー,'
        'リンク長,リンク種別,自専道,通行可不可,幅員区分,車線数,交通規制,'
        '現・旧区分,基本区間,左右フラグ,ﾕﾆｰｸｺｰﾄﾞ,所要時間（0.01分）'
    )

    link_cnt     = 0
    feat_cnt     = 0
    skip_loop    = 0
    link_id      = 0
    parquet_rows = []

    with OUT_LINKS_CSV.open('w', encoding='utf-8-sig') as csv_out:
        csv_out.write(CSV_HEADER + '\n')

        for fi, fp in enumerate(geojson_files, 1):
            with fp.open(encoding='utf-8') as f:
                data = json.load(f)
            features = data['features']

            for feat in features:
                if USE_FILTER and not passes_filter(feat['properties']):
                    continue
                feat_cnt += 1

                props  = feat['properties']
                coords = feat['geometry']['coordinates']

                n13_002 = str(props.get('N13_002', '1'))
                n13_003 = str(props.get('N13_003', '3'))
                n13_004 = str(props.get('N13_004', '0'))
                n13_005 = str(props.get('N13_005', '0'))
                n13_006 = str(props.get('N13_006', '0'))
                n13_007 = str(props.get('N13_007', '0'))
                n13_008 = str(props.get('N13_008', '0'))

                pt1 = (round(coords[0][0],  9), round(coords[0][1],  9))
                pt2 = (round(coords[-1][0], 9), round(coords[-1][1], 9))

                n1 = coord_to_nid[pt1]
                n2 = coord_to_nid[pt2]

                if n1 == n2:
                    skip_loop += 1
                    continue

                dist_m = max(1, round(line_length_m(coords)))
                speed  = WALK_KMH if MODE == 'walk' else SPEED_KMH.get(n13_003, 20)
                tim    = max(1, round(dist_m / speed * 6.0))
                rtype  = ROAD_TYPE.get(n13_003, 7)
                jis    = JISEN.get(n13_003, 2)
                mesh   = mesh2_code(pt1[0], pt1[1])

                link_id += 1
                uid_f = f'{link_id:016d}'
                uid_r = f'R{link_id:015d}'

                common_cols = (f'0,{rtype},{n13_008},0,{dist_m},'
                               f'{rtype},{jis},0,{n13_006},{n13_007},0,0,00000000000')

                # CSV: 正方向
                csv_out.write(f'{mesh},{n1},{n2},{common_cols},1,{uid_f},{tim}\n')
                link_cnt += 1

                # CSV: 逆方向（全道路・KSJに一方通行フィールドなし）
                csv_out.write(f'{mesh},{n2},{n1},{common_cols},2,{uid_r},{tim}\n')
                link_cnt += 1

                # Parquet用: 正方向のみ（全補間点を含む座標列）
                parquet_rows.append({
                    'node1'       : n1,
                    'node2'       : n2,
                    'mesh2'       : mesh,
                    'N13_002'     : n13_002,
                    'N13_003'     : n13_003,
                    'N13_004'     : n13_004,
                    'N13_005'     : n13_005,
                    'N13_006'     : n13_006,
                    'N13_007'     : n13_007,
                    'N13_008'     : n13_008,
                    'dist_m'      : dist_m,
                    'time_001min' : tim,
                    'road_type'   : rtype,
                    'geometry'    : LineString(coords),
                })

                if feat_cnt % PROG_STEP == 0:
                    pct = feat_cnt * 100 // max(pass1_filtered, 1)
                    print(f'  {feat_cnt:>8,} / {pass1_filtered:,} ({pct}%)  '
                          f'リンク: {link_cnt:,}  ({time.time()-t0:.0f}s)')

    print(f'  出力リンク数(CSV): {link_cnt:,}')
    print(f'  自己ループ除去   : {skip_loop:,}')
    print(f'  経過時間: {time.time()-t0:.0f}s')

    # ── [3/4] リンクParquet出力 ───────────────────────────────
    print(f'\n[3/4] リンクParquet出力 ...')
    gdf_links = gpd.GeoDataFrame(parquet_rows, crs=CRS)
    del parquet_rows
    gdf_links.to_parquet(OUT_LINKS_PARQUET)
    print(f'      完了  ({time.time()-t0:.0f}s)')

    # ── [4/4] ノードCSV + ノードParquet 出力 ────────────────
    print(f'\n[4/4] ノードCSV・ノードParquet出力 ...')
    items      = list(coord_to_nid.items())
    node_total = len(items)

    with OUT_NODES_CSV.open('w', encoding='utf-8-sig') as ncsv:
        ncsv.write('node_id,lon,lat\n')
        for (lon, lat), node_id in items:
            ncsv.write(f'{node_id},{lon},{lat}\n')

    node_rows = [
        {'node_id': nid, 'geometry': Point(lon, lat)}
        for (lon, lat), nid in items
    ]
    gdf_nodes = gpd.GeoDataFrame(node_rows, crs=CRS)
    del node_rows
    gdf_nodes.to_parquet(OUT_NODES_PARQUET)
    print(f'      完了  ({time.time()-t0:.0f}s)')

    # ── [5/4] クリップ（--city / --pref）+ QML出力 ────────────
    print(f'\n[5/4] クリップ・QML出力 ...')
    OUT_LINKS_QML = OUT_DIR / f'KSJ_N13-24_{CASE_NAME}_道路リンク.qml'

    clip_gdf = None
    if CLIP_CITY:
        city_names = [c.strip() for c in CLIP_CITY.split(',')]
        print(f'      市区町村クリップ: {", ".join(city_names)}')
        raw_gdf = gpd.read_parquet(CITY_PARQUET)
        if CLIP_PREF:
            pref_names = [p.strip() for p in CLIP_PREF.split(',')]
            raw_gdf = raw_gdf[raw_gdf['N03_001'].isin(pref_names)]
        clip_gdf = raw_gdf[raw_gdf['N03_004'].isin(city_names)]
        missing = [c for c in city_names if c not in raw_gdf['N03_004'].values]
        if missing:
            print(f'      警告: 見つからない市区町村名: {missing}')
    elif CLIP_PREF:
        pref_names = [p.strip() for p in CLIP_PREF.split(',')]
        print(f'      都道府県クリップ: {", ".join(pref_names)}')
        pref_df = gpd.read_parquet(PREF_PARQUET)
        clip_gdf = pref_df[pref_df['prefecture'].isin(pref_names)]
        missing = [p for p in pref_names if p not in pref_df['prefecture'].values]
        if missing:
            print(f'      警告: 見つからない都道府県名: {missing}')

    if clip_gdf is not None and len(clip_gdf) > 0:
        import shapely as shp
        from shapely.ops import unary_union
        clip_poly = unary_union(clip_gdf.to_crs(CRS).geometry)

        # sindex でバウンディングボックス事前フィルター → shapely ufunc で正確フィルター（高速）
        n_links_orig = len(gdf_links)
        cands = list(gdf_links.sindex.intersection(clip_poly.bounds))
        if cands:
            sub = gdf_links.iloc[cands]
            mask = shp.intersects(sub.geometry.values, clip_poly)
            gdf_links = sub.iloc[mask]
        else:
            gdf_links = gdf_links.iloc[[]]
        print(f'      リンクParquetクリップ: {len(gdf_links):,} / {n_links_orig:,}')
        gdf_links.to_parquet(OUT_LINKS_PARQUET)

        n_nodes_orig = len(gdf_nodes)
        cands = list(gdf_nodes.sindex.intersection(clip_poly.bounds))
        if cands:
            sub = gdf_nodes.iloc[cands]
            mask = shp.within(sub.geometry.values, clip_poly)
            gdf_nodes = sub.iloc[mask]
        else:
            gdf_nodes = gdf_nodes.iloc[[]]
        print(f'      ノードParquetクリップ: {len(gdf_nodes):,} / {n_nodes_orig:,}')
        gdf_nodes.to_parquet(OUT_NODES_PARQUET)

    write_road_links_qml(OUT_LINKS_QML)
    print(f'      QML: {OUT_LINKS_QML.name}')
    del gdf_links, gdf_nodes
    print(f'      完了  ({time.time()-t0:.0f}s)')

    # ── 完了サマリー ──────────────────────────────────────────
    elapsed = time.time() - t0
    print(f'\n===== 完了 ({elapsed:.0f}s) =====')
    print(f'  対象範囲                 : {mesh_desc}')
    print(f'  フィルター後フィーチャ数 : {feat_cnt:,}')
    print(f'  ユニークノード数         : {node_total:,}')
    print(f'  使用2次メッシュ数        : {len(mesh_counters):,}')
    print(f'  出力リンク数(CSV)        : {link_cnt:,}')
    print(f'  自己ループ除去           : {skip_loop:,}')
    print()
    for p in [OUT_LINKS_CSV, OUT_LINKS_PARQUET, OUT_LINKS_QML, OUT_NODES_CSV, OUT_NODES_PARQUET]:
        if p.exists():
            size_mb = os.path.getsize(p) / 1024 / 1024
            print(f'  {p.name:<50s}  {size_mb:6.0f} MB')


if __name__ == '__main__':
    main()
