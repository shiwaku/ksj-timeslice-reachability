#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KSJ道路ネットワーク アクセスリンク生成スクリプト

入力 :
  KSJ_N13-24_XXXX_道路ノード.csv  ノードCSV（node_id, lon, lat）
出力（XXXX=1次メッシュコードまたはケース名, N=レベル番号）:
  KSJ_N13-24_XXXX_アクセスリンク_LN.csv          アクセスリンクCSV（DRM3003形式準拠）
  KSJ_N13-24_XXXX_アクセスリンク_make_net_LN.csv  make_net.f用3列CSV
  KSJ_N13-24_XXXX_アクセスリンク_LN.geojson      アクセスリンクGeoJSON（可視化用）
  KSJ_N13-24_XXXX_アクセスリンク_LN.parquet      アクセスリンクGeoParquet（QGIS・分析用）

使い方:
  # 単一1次メッシュ
  python3 make_access_links.py --mesh 5238 --case 5238_shizuoka --level 4          # L4・車
  python3 make_access_links.py --mesh 5238 --case 5238_shizuoka_walk --level 6 --mode walk  # L6・徒歩

  # 複数1次メッシュ（地方ブロック / 都道府県）・海上除外自動適用
  python3 make_access_links.py \\
    --meshes 5239,5240,5338,5339,5340,5438,5439,5440,5539,5540 \\
    --case kanto_block --level 4 \\
    --pref 茨城県,栃木県,群馬県,埼玉県,千葉県,東京都,神奈川県,山梨県  # 関東ブロック（L4・都県境クリップ）
  python3 make_access_links.py \\
    --meshes 5137,5138,5237,5238,5337,5338 --case shizuoka_pref --level 5  # 都道府県（L5）

  # 全国（L3・1kmメッシュ）
  python3 make_access_links.py --nationwide --level 3 --case nationwide

  # 都道府県境界でクリップ（--pref）
  python3 make_access_links.py \\
    --meshes 5137,5138,5237,5238,5337,5338 --case shizuoka_pref_clip --level 5 --pref 静岡県

  # 任意ポリゴンでクリップ（--clip-geom）
  python3 make_access_links.py \\
    --mesh 5238 --case shizuoka_city --level 6 --clip-geom /path/to/shizuoka_city.geojson

  # census SHPモード（100mメッシュ重心を使用）:
  python3 make_access_links.py \
    --mesh 5238 --case 5238_yaizu_walk --mode walk \
    --census-shp /path/to/census.shp \
    --extra-points /path/to/shelters.geojson --extra-filter "津波=1" \
    --nodes-csv /path/to/道路ノード.csv

【速度モード（--mode）】
  vehicle（デフォルト）: 5km/h（ゆっくりした歩行想定）
  walk: 3.6km/h（津波避難徒歩速度）

【アクセスリンクとは】
  各メッシュ重心（仮想ノード）→ 最近傍道路ノードをつなぐリンク。
  経路探索の出発・到着点として使用する。

【メッシュ仕様】
  レベル3（1km） :  8桁コード, 1次メッシュあたり   6,400個
  レベル4（500m）:  9桁コード, 1次メッシュあたり  25,600個
  レベル5（250m）: 10桁コード, 1次メッシュあたり 102,400個
  レベル6（125m）: 11桁コード, 1次メッシュあたり 409,600個
  census 100m    : 10桁コード（統計メッシュ方式）, SHPから読み込み

【アクセスノードID】
  道路ノードID（先頭2桁がPP=30-68）との衝突を避けるため先頭で区別:
    レベル3: '000' + メッシュコード(8桁)  → '000XXXXXXXX'（11桁）
    レベル4: '00'  + メッシュコード(9桁)  → '00XXXXXXXXX'（11桁）
    レベル5: '0'   + メッシュコード(10桁) → '0XXXXXXXXXX'（11桁）
    レベル6: '0'   + 10桁suffix          → '0XXXXXXXXXX'（11桁）
    census : '0'   + 9桁suffix           → '0XXXXXXXXX' （10桁, I11で2桁パディング）

【--meshes モード（複数1次メッシュ）】
  --meshes 5338,5339,5340,... と指定すると複数1次メッシュを横断してメッシュを列挙する。
  prefecture.parquet（02_ShortestRouteSearch/out/prefecture.parquet）を使って
  海上メッシュ（都道府県ポリゴン外）を自動除外する。
  --nodes-csv を省略した場合、KSJ_N13-24_{CASE_NAME}_道路ノード.csv を自動参照する。

【リンク仕様】
  方向: 双方向（往復2本）
  速度: ACCESS_KMH（デフォルト5km/h = 徒歩）
  リンク長: 重心→最近傍道路ノード間のHaversine距離
"""

import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.spatial import KDTree

# ============================
# 引数解析
# ============================
parser = argparse.ArgumentParser(description='KSJ道路ネットワーク アクセスリンク生成（統合版）')
parser.add_argument('--mesh',  default=None,
                    help='単一1次メッシュコード（例: 5238）')
parser.add_argument('--meshes', default=None,
                    help='複数1次メッシュコード（カンマ区切り、例: 5338,5339,5340）。--case 必須')
parser.add_argument('--nationwide', action='store_true',
                    help='全国全メッシュ（L3専用）。--level 3 と組み合わせて使用')
parser.add_argument('--case',  default=None,
                    help='ケース名（出力サブディレクトリ名）。省略時: --mesh → メッシュコード, --nationwide → nationwide')
parser.add_argument('--level', type=int, default=4, choices=[3, 4, 5, 6],
                    help='メッシュレベル: 3=1km, 4=500m（デフォルト）, 5=250m, 6=125m。--census-shp指定時は無視')
parser.add_argument('--mode',  default='vehicle', choices=['vehicle', 'walk'],
                    help='速度モード: vehicle=5km/h（デフォルト）, walk=3.6km/h（徒歩速度）')
parser.add_argument('--walk-kmh', type=float, default=None,
                    help='walk モード時の速度（km/h）。省略時は mode=walk なら3.6、vehicle なら5.0')
# ── 空間フィルター ──
parser.add_argument('--pref', default=None,
                    help='都道府県名でメッシュをクリップ（カンマ区切りで複数指定可。prefecture.parquetを使用）'
                         '例: 静岡県  / 茨城県,栃木県,群馬県,埼玉県,千葉県,東京都,神奈川県,山梨県')
parser.add_argument('--clip-geom', default=None,
                    help='任意ポリゴンファイル（GeoJSON/GeoParquet）でメッシュをクリップ')
parser.add_argument('--city', default=None,
                    help='市区町村名でメッシュをクリップ（カンマ区切りで複数指定可。N03-20250101.geojsonを使用）'
                         '例: 高松市 / 吉岡町（--prefと組み合わせで都道府県絞り込み可）')
# --- census SHP モード ---
parser.add_argument('--census-shp', default=None,
                    help='国勢調査100mメッシュSHPファイル。指定時はSHP重心モード（L6列挙の代わりに使用）')
parser.add_argument('--extra-points', default=None,
                    help='追加ポイントGeoJSON（--census-shpと組み合わせ、避難場所等を追加）')
parser.add_argument('--extra-filter', default=None,
                    help='追加ポイントフィルター: "field=value" 形式 (例: "津波=1")')
parser.add_argument('--nodes-csv', default=None,
                    help='道路ノードCSVパス（省略時は --case ディレクトリ内を自動検索）')
parser.add_argument('--out-dir', default=None,
                    help='出力先ディレクトリ（省略時: 01_MakeNetwork/{case}/）')
args = parser.parse_args()

# ── 対象範囲モードの決定 ──
NATIONWIDE_MODE = args.nationwide
MULTI_MESH_MODE = args.meshes is not None

if NATIONWIDE_MODE:
    MESH_LIST = None   # 動的取得（nationwide）
    MESH_CODE = None
    if args.level != 3:
        print('警告: --nationwide は L3（--level 3）専用です。--level 3 を自動設定します。')
        args.level = 3
elif MULTI_MESH_MODE:
    MESH_LIST = [m.strip() for m in args.meshes.split(',')]
    MESH_CODE = MESH_LIST[0]
    if not args.case:
        parser.error('--meshes 使用時は --case が必須です')
else:
    # 単一メッシュ（--mesh が未指定の場合はエラー）
    if not args.mesh:
        parser.error('--mesh, --meshes, --nationwide のいずれかが必要です')
    MESH_LIST = [args.mesh]
    MESH_CODE = args.mesh

# ── ケース名デフォルト ──
if args.case:
    CASE_NAME = args.case
elif NATIONWIDE_MODE:
    CASE_NAME = 'nationwide'
elif args.mesh:
    CASE_NAME = args.mesh
else:
    CASE_NAME = MESH_LIST[0]   # fallback

LEVEL       = args.level
MODE        = args.mode
CENSUS_MODE = args.census_shp is not None
PREF_NAME   = args.pref
CLIP_GEOM   = args.clip_geom
CITY_NAME   = args.city

# ============================
# 設定
# ============================
BASE_DIR   = Path(__file__).parent.parent
OUT_DIR    = Path(args.out_dir) if args.out_dir else BASE_DIR / 'network' / CASE_NAME
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 出力ファイル名プレフィクス（常にケース名を使用）
_prefix = CASE_NAME

# 道路ノードCSV（--nodes-csv 優先、なければケース名から自動解決）
if args.nodes_csv:
    IN_NODES_CSV = Path(args.nodes_csv)
else:
    IN_NODES_CSV = OUT_DIR / f'KSJ_N13-24_{CASE_NAME}_道路ノード.csv'

# 出力ファイル名
_suffix = 'census' if CENSUS_MODE else f'L{LEVEL}'
OUT_ACCESS_CSV     = OUT_DIR / f'KSJ_N13-24_{_prefix}_アクセスリンク_{_suffix}.csv'
OUT_ACCESS_MAKENET = OUT_DIR / f'KSJ_N13-24_{_prefix}_アクセスリンク_make_net_{_suffix}.csv'
OUT_ACCESS_GEOJSON = OUT_DIR / f'KSJ_N13-24_{_prefix}_アクセスリンク_{_suffix}.geojson'

if args.walk_kmh is not None:
    ACCESS_KMH = args.walk_kmh
else:
    ACCESS_KMH = 3.6 if MODE == 'walk' else 5    # アクセスリンク速度（km/h）

# prefecture.parquet のパス（複数メッシュ海上除外用）
PREF_PARQUET = BASE_DIR / 'data' / 'prefecture.parquet'
# 市区町村ポリゴン（--city クリップ用）
CITY_PARQUET = BASE_DIR / 'data' / 'city.parquet'

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


def mesh4_code(lon, lat):
    """経緯度 → 4次メッシュコード（9桁）"""
    p = int(lat * 1.5)
    u = int(lon - 100)
    q = int((lat - p / 1.5) * 12)
    r = int((lon - (100 + u)) * 8)
    lat_rem2 = (lat - p / 1.5) - q / 12
    lon_rem2 = (lon - (100 + u)) - r / 8
    s = int(lat_rem2 * 120)
    t = int(lon_rem2 * 80)
    lat_rem3 = lat_rem2 - s / 120
    lon_rem3 = lon_rem2 - t / 80
    v_row = int(lat_rem3 * 240)
    v_col = int(lon_rem3 * 160)
    v = v_row * 2 + v_col + 1
    return f'{p:02d}{u:02d}{q}{r}{s}{t}{v}'


def mesh4_centroid(code):
    """4次メッシュコード(9桁) → 重心(lon, lat)"""
    p = int(code[0:2]); u = int(code[2:4])
    q = int(code[4]);   r = int(code[5])
    s = int(code[6]);   t = int(code[7])
    v = int(code[8])
    lat_base = p / 1.5 + q / 12 + s / 120
    lon_base = 100 + u + r / 8 + t / 80
    d_lat = 1 / 240
    d_lon = 1 / 160
    v_row = (v - 1) // 2
    v_col = (v - 1) % 2
    lat = lat_base + v_row * d_lat + d_lat / 2
    lon = lon_base + v_col * d_lon + d_lon / 2
    return lon, lat


def mesh5_centroid(code):
    """5次メッシュコード(10桁) → 重心(lon, lat)"""
    p = int(code[0:2]); u = int(code[2:4])
    q = int(code[4]);   r = int(code[5])
    s = int(code[6]);   t = int(code[7])
    v = int(code[8]);   w = int(code[9])
    lat_base = p / 1.5 + q / 12 + s / 120
    lon_base = 100 + u + r / 8 + t / 80
    d_lat4 = 1 / 240;  d_lon4 = 1 / 160
    d_lat5 = d_lat4 / 2; d_lon5 = d_lon4 / 2
    v_row = (v - 1) // 2; v_col = (v - 1) % 2
    w_row = (w - 1) // 2; w_col = (w - 1) % 2
    lat = lat_base + v_row * d_lat4 + w_row * d_lat5 + d_lat5 / 2
    lon = lon_base + v_col * d_lon4 + w_col * d_lon5 + d_lon5 / 2
    return lon, lat


def mesh6_centroid(code):
    """6次メッシュコード(11桁) → 重心(lon, lat)"""
    p = int(code[0:2]); u = int(code[2:4])
    q = int(code[4]);   r = int(code[5])
    s = int(code[6]);   t = int(code[7])
    v = int(code[8]);   w = int(code[9]); x = int(code[10])
    lat_base = p / 1.5 + q / 12 + s / 120
    lon_base = 100 + u + r / 8 + t / 80
    d_lat4 = 1 / 240;  d_lon4 = 1 / 160
    d_lat5 = d_lat4 / 2; d_lon5 = d_lon4 / 2
    d_lat6 = d_lat5 / 2; d_lon6 = d_lon5 / 2
    v_row = (v - 1) // 2; v_col = (v - 1) % 2
    w_row = (w - 1) // 2; w_col = (w - 1) % 2
    x_row = (x - 1) // 2; x_col = (x - 1) % 2
    lat = lat_base + v_row * d_lat4 + w_row * d_lat5 + x_row * d_lat6 + d_lat6 / 2
    lon = lon_base + v_col * d_lon4 + w_col * d_lon5 + x_col * d_lon6 + d_lon6 / 2
    return lon, lat


def mesh100m_code(lon, lat):
    """経緯度 → 100mメッシュコード（10桁・統計メッシュ方式・国勢調査MESH_CODEと同形式）

    3次メッシュ（1km）を緯度10分割 × 経度10分割した100mメッシュの識別コード。
    L6（125m・2進分割）とは異なるグリッド体系。
    """
    p = int(lat * 1.5)
    u = int(lon - 100)
    q = int((lat - p / 1.5) * 12)
    r = int((lon - (100 + u)) * 8)
    lat_rem2 = (lat - p / 1.5) - q / 12
    lon_rem2 = (lon - (100 + u)) - r / 8
    s = int(lat_rem2 * 120)
    t = int(lon_rem2 * 80)
    lat_rem3 = lat_rem2 - s / 120
    lon_rem3 = lon_rem2 - t / 80
    x = int(lat_rem3 * 1200)   # 0-9: 3次メッシュ内の緯度方向位置
    y = int(lon_rem3 * 800)    # 0-9: 3次メッシュ内の経度方向位置
    return f'{p:02d}{u:02d}{q}{r}{s}{t}{x}{y}'


def enum_mesh4_in_mesh1(mesh1_code):
    """1次メッシュコード(4桁) → 内包する全4次メッシュコードを列挙"""
    p = int(mesh1_code[0:2])
    u = int(mesh1_code[2:4])
    codes = []
    for q in range(8):
        for r in range(8):
            for s in range(10):
                for t in range(10):
                    for v in range(1, 5):
                        codes.append(f'{p:02d}{u:02d}{q}{r}{s}{t}{v}')
    return codes


def enum_mesh5_in_mesh1(mesh1_code):
    """1次メッシュコード(4桁) → 内包する全5次メッシュコードを列挙"""
    p = int(mesh1_code[0:2])
    u = int(mesh1_code[2:4])
    codes = []
    for q in range(8):
        for r in range(8):
            for s in range(10):
                for t in range(10):
                    for v in range(1, 5):
                        for w in range(1, 5):
                            codes.append(f'{p:02d}{u:02d}{q}{r}{s}{t}{v}{w}')
    return codes


def enum_mesh6_in_mesh1(mesh1_code):
    """1次メッシュコード(4桁) → 内包する全6次メッシュコードを列挙"""
    p = int(mesh1_code[0:2])
    u = int(mesh1_code[2:4])
    codes = []
    for q in range(8):
        for r in range(8):
            for s in range(10):
                for t in range(10):
                    for v in range(1, 5):
                        for w in range(1, 5):
                            for x in range(1, 5):
                                codes.append(f'{p:02d}{u:02d}{q}{r}{s}{t}{v}{w}{x}')
    return codes


def filter_land_meshes(codes, centroid_fn):
    """メッシュコードリストから海上メッシュを除外して陸上メッシュのみ返す。

    prefecture.parquet（都道府県境界ポリゴン）との空間結合で判定する。
    geopandas の sjoin を使ったベクトル化処理で高速に処理する。

    Args:
        codes:        メッシュコードリスト
        centroid_fn:  メッシュコード → (lon, lat) 関数

    Returns:
        land_codes:   陸上メッシュのみのリスト（元の順序を保持）
        n_sea:        除外された海上メッシュ数
    """
    import geopandas as gpd
    import pandas as pd
    from shapely.geometry import Point

    print(f'      海上除外: prefecture.parquet 読み込み中 ...')
    if not PREF_PARQUET.exists():
        print(f'      警告: {PREF_PARQUET} が見つかりません。海上除外をスキップします。')
        return codes, 0

    pref_gdf = gpd.read_parquet(PREF_PARQUET)
    # CRS を EPSG:4326 に統一（EPSG:6668 ≈ EPSG:4326 だが明示的に揃える）
    if pref_gdf.crs is None:
        pref_gdf = pref_gdf.set_crs('EPSG:6668')
    pref_gdf = pref_gdf.to_crs('EPSG:4326')

    print(f'      海上除外: {len(codes):,} メッシュの重心計算中 ...')
    lons = []
    lats = []
    for code in codes:
        lon, lat = centroid_fn(code)
        lons.append(lon)
        lats.append(lat)

    centroid_gdf = gpd.GeoDataFrame(
        {'code': codes},
        geometry=gpd.points_from_xy(lons, lats),
        crs='EPSG:4326'
    )

    print(f'      海上除外: sjoin で陸上判定中 ...')
    joined = centroid_gdf.sjoin(
        pref_gdf[['geometry']],
        how='left',
        predicate='within'
    )
    # sjoin は1対多になる場合があるため重複排除（元インデックスで1件に絞る）
    land_idx = set(joined[joined['index_right'].notna()].index.unique())

    land_codes = [code for i, code in enumerate(codes) if i in land_idx]
    n_sea      = len(codes) - len(land_codes)
    print(f'      海上除外: 陸上 {len(land_codes):,} / 全体 {len(codes):,}（除外 {n_sea:,} 件）')
    return land_codes, n_sea


def filter_by_polygon(codes, centroid_fn, poly_gdf):
    """メッシュコードリストを任意ポリゴン（GeoDataFrame）でクリップして返す。

    メッシュ重心が poly_gdf のいずれかのポリゴン内に含まれるものだけを残す。
    filter_land_meshes() の後に適用することで行政界によるクリップが可能。

    Args:
        codes:       メッシュコードリスト
        centroid_fn: メッシュコード → (lon, lat) 関数
        poly_gdf:    クリップ用 GeoDataFrame（EPSG:4326 に変換済みを期待）

    Returns:
        clipped_codes: クリップ後メッシュコードリスト
        n_removed:     除外されたメッシュ数
    """
    import geopandas as gpd

    print(f'      ポリゴンクリップ: {len(codes):,} メッシュの重心計算中 ...')
    lons, lats = [], []
    for code in codes:
        lon, lat = centroid_fn(code)
        lons.append(lon)
        lats.append(lat)

    centroid_gdf = gpd.GeoDataFrame(
        {'code': codes},
        geometry=gpd.points_from_xy(lons, lats),
        crs='EPSG:4326'
    )
    if poly_gdf.crs is None:
        poly_gdf = poly_gdf.set_crs('EPSG:4326')
    else:
        poly_gdf = poly_gdf.to_crs('EPSG:4326')

    joined = centroid_gdf.sjoin(poly_gdf[['geometry']], how='left', predicate='within')
    in_idx = set(joined[joined['index_right'].notna()].index.unique())

    clipped = [code for i, code in enumerate(codes) if i in in_idx]
    n_removed = len(codes) - len(clipped)
    print(f'      ポリゴンクリップ: {len(clipped):,} / {len(codes):,}（除外 {n_removed:,} 件）')
    return clipped, n_removed


def mesh3_centroid(code):
    """3次メッシュコード（8桁）→ 重心 (lon, lat)"""
    p = int(code[0:2]); u = int(code[2:4])
    q = int(code[4]);   r = int(code[5])
    s = int(code[6]);   t = int(code[7])
    lat_base = p / 1.5 + q / 12 + s / 120
    lon_base = 100 + u + r / 8 + t / 80
    d_lat3 = 1 / 120   # 3次メッシュの緯度幅（≒925m）
    d_lon3 = 1 / 80    # 3次メッシュの経度幅（≒1039m at 35°N）
    return lon_base + d_lon3 / 2, lat_base + d_lat3 / 2


def enum_mesh3_in_mesh1(mesh1_code):
    """1次メッシュコード（4桁）→ 内包する全3次メッシュコード（8桁）を列挙"""
    p = int(mesh1_code[0:2])
    u = int(mesh1_code[2:4])
    codes = []
    for q in range(8):
        for r in range(8):
            for s in range(10):
                for t in range(10):
                    codes.append(f'{p:02d}{u:02d}{q}{r}{s}{t}')
    return codes


# ============================
# メイン処理
# ============================
def main():
    t0 = time.time()

    # ── [1/4] 道路ノードCSV 読み込み ────────────────────────
    print('[1/4] 道路ノードCSV読み込み中 ...')
    print(f'      {IN_NODES_CSV}')
    node_ids = []
    lons = []
    lats = []
    with IN_NODES_CSV.open(encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            node_ids.append(row['node_id'])
            lons.append(float(row['lon']))
            lats.append(float(row['lat']))
    print(f'      道路ノード数: {len(node_ids):,}  ({time.time()-t0:.1f}s)')

    # ── [2/4] KDTree 構築 ────────────────────────────────────
    print('\n[2/4] KDTree構築中 ...')
    # 緯度経度をラジアンに変換して球面KDTree（近似：平面投影）
    # 日本国内の小距離なら平面近似で十分
    coords = np.column_stack([lons, lats])
    tree = KDTree(coords)
    print(f'      完了  ({time.time()-t0:.1f}s)')

    # ── [3/4] メッシュ列挙 & 最近傍探索 ─────────────────────
    print(f'\n[3/4] アクセスリンク生成中 ...')
    print(f'      CSV         : {OUT_ACCESS_CSV.name}')
    print(f'      make_net用  : {OUT_ACCESS_MAKENET.name}')
    print(f'      GeoJSON     : {OUT_ACCESS_GEOJSON.name}')

    if CENSUS_MODE:
        # ---- census SHP モード: 100mメッシュ重心を使用 ----
        import geopandas as gpd
        import warnings
        import json as _json

        print(f'      モード: census SHP（100mメッシュ重心）')
        census_df = gpd.read_file(args.census_shp, encoding='cp932').to_crs(epsg=4326)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            census_df = census_df.copy()
            census_df['_cx'] = census_df.geometry.centroid.x
            census_df['_cy'] = census_df.geometry.centroid.y

        mesh_entries = []  # (c_lon, c_lat, mesh_code_str, access_nid_str)
        seen_nids = set()
        for _, row in census_df.iterrows():
            mesh_code_str = str(row['MESH_CODE'])
            access_nid = f'0{mesh_code_str[1:]}'   # "0" + 9桁suffix = 10桁
            if access_nid in seen_nids:
                continue
            seen_nids.add(access_nid)
            mesh_entries.append((float(row['_cx']), float(row['_cy']), mesh_code_str, access_nid))
        print(f'      SHPメッシュ数: {len(mesh_entries):,}')

        # 追加ポイント（避難場所等）
        if args.extra_points:
            ef, ev = args.extra_filter.split('=', 1) if args.extra_filter else (None, None)
            with open(args.extra_points, encoding='utf-8') as fp:
                gj_extra = _json.load(fp)
            extra_added = 0
            for feat in gj_extra['features']:
                if ef and str(feat['properties'].get(ef)) != ev:
                    continue
                lon, lat = feat['geometry']['coordinates'][:2]
                code = mesh100m_code(lon, lat)
                nid = f'0{code[1:]}'
                if nid in seen_nids:
                    continue
                seen_nids.add(nid)
                # 実座標を使用（最近傍道路ノード探索の精度向上）
                mesh_entries.append((lon, lat, code, nid))
                extra_added += 1
            print(f'      追加ポイント（新規）: {extra_added:,}件')

        total = len(mesh_entries)
        print(f'      合計: {total:,}エントリ')

    else:
        # ---- L3/L4/L5/L6 列挙モード ----
        mesh_size_str = {3: '1km', 4: '500m', 5: '250m', 6: '125m'}[LEVEL]
        enum_fn     = {3: enum_mesh3_in_mesh1,
                       4: enum_mesh4_in_mesh1,
                       5: enum_mesh5_in_mesh1,
                       6: enum_mesh6_in_mesh1}[LEVEL]
        centroid_fn = {3: mesh3_centroid,
                       4: mesh4_centroid,
                       5: mesh5_centroid,
                       6: mesh6_centroid}[LEVEL]

        # 対象1次メッシュリストの決定
        if NATIONWIDE_MODE:
            # 全国: input/ から1次メッシュコードを動的取得
            geojson_dir = Path(__file__).parent.parent / 'input'
            geojson_files = sorted(geojson_dir.glob('*/N13-24_*.geojson'))
            target_mesh1 = [
                fp.stem.replace('N13-24_', '')
                for fp in geojson_files
            ]
            print(f'      モード: 全国・L{LEVEL}（{mesh_size_str}メッシュ）・{len(target_mesh1)}1次メッシュ')
        else:
            target_mesh1 = MESH_LIST
            range_desc = '複数1次メッシュ' if MULTI_MESH_MODE else '単一1次メッシュ'
            print(f'      モード: L{LEVEL}（{mesh_size_str}メッシュ）・{range_desc}: {target_mesh1}')

        raw_codes = []
        for m1 in target_mesh1:
            raw_codes.extend(enum_fn(m1))
        print(f'      列挙メッシュ数（フィルター前）: {len(raw_codes):,}')

        # 海上除外（--nationwide または --meshes の場合）
        if NATIONWIDE_MODE or MULTI_MESH_MODE:
            raw_codes, _ = filter_land_meshes(raw_codes, centroid_fn)

        # 都道府県クリップ（--pref, カンマ区切りで複数可）
        if PREF_NAME and not CITY_NAME:
            import geopandas as gpd
            pref_names = [p.strip() for p in PREF_NAME.split(',')]
            print(f'      都道府県クリップ: {", ".join(pref_names)}')
            pref_gdf = gpd.read_parquet(PREF_PARQUET)
            pref_gdf = pref_gdf[pref_gdf['prefecture'].isin(pref_names)]
            missing = [p for p in pref_names if p not in pref_gdf['prefecture'].values]
            if missing:
                print(f'      警告: 見つからない都道府県名: {missing}')
            if len(pref_gdf) == 0:
                print(f'      警告: 一致する都道府県が見つかりません。クリップをスキップします。')
            else:
                raw_codes, _ = filter_by_polygon(raw_codes, centroid_fn, pref_gdf)

        # 市区町村クリップ（--city）
        if CITY_NAME:
            import geopandas as gpd
            city_names = [c.strip() for c in CITY_NAME.split(',')]
            print(f'      市区町村クリップ: {", ".join(city_names)}')
            city_gdf = gpd.read_parquet(CITY_PARQUET)
            if PREF_NAME:
                pref_names = [p.strip() for p in PREF_NAME.split(',')]
                city_gdf = city_gdf[city_gdf['N03_001'].isin(pref_names)]
            city_gdf = city_gdf[city_gdf['N03_004'].isin(city_names)]
            missing = [c for c in city_names if c not in city_gdf['N03_004'].values]
            if missing:
                print(f'      警告: 見つからない市区町村名: {missing}')
            if len(city_gdf) == 0:
                print(f'      警告: 一致する市区町村が見つかりません。クリップをスキップします。')
            else:
                city_gdf = city_gdf.to_crs('EPSG:4326')
                raw_codes, _ = filter_by_polygon(raw_codes, centroid_fn, city_gdf)

        # 任意ポリゴンクリップ（--clip-geom）
        if CLIP_GEOM:
            import geopandas as gpd
            print(f'      ポリゴンクリップ: {CLIP_GEOM}')
            clip_path = Path(CLIP_GEOM)
            if clip_path.suffix == '.parquet':
                clip_gdf = gpd.read_parquet(clip_path)
            else:
                clip_gdf = gpd.read_file(clip_path)
            raw_codes, _ = filter_by_polygon(raw_codes, centroid_fn, clip_gdf)

        total = len(raw_codes)
        print(f'      {LEVEL}次メッシュ数（最終）: {total:,}')

    CSV_HEADER = (
        '2次ﾒｯｼｭ,ﾉｰﾄﾞ1,ﾉｰﾄﾞ2,管理者,道路種別,路線番号,市町村コー,'
        'リンク長,リンク種別,自専道,通行可不可,幅員区分,車線数,交通規制,'
        '現・旧区分,基本区間,左右フラグ,ﾕﾆｰｸｺｰﾄﾞ,所要時間（0.01分）'
    )

    link_cnt = 0
    first_feature = True

    with (OUT_ACCESS_CSV.open('w', encoding='utf-8-sig') as csv_out,
          OUT_ACCESS_MAKENET.open('w', encoding='utf-8-sig') as mn_out,
          OUT_ACCESS_GEOJSON.open('w', encoding='utf-8') as gjson_out):

        csv_out.write(CSV_HEADER + '\n')
        mn_out.write(f'アクセスノードID,道路ノードID,距離\n')
        gjson_out.write('{"type":"FeatureCollection","features":[\n')

        if CENSUS_MODE:
            iterator = enumerate(mesh_entries)
        else:
            iterator = enumerate(raw_codes)

        for i, item in iterator:
            if CENSUS_MODE:
                c_lon, c_lat, mcode, access_nid = item
            else:
                mcode = item
                c_lon, c_lat = centroid_fn(mcode)
                # アクセスノードID（道路ノードIDと非衝突）
                if LEVEL == 3:
                    access_nid = f'000{mcode}'     # "000" + 8桁 = 11桁
                elif LEVEL == 4:
                    access_nid = f'00{mcode}'      # "00"  + 9桁 = 11桁
                elif LEVEL == 5:
                    access_nid = f'0{mcode}'       # "0"   + 10桁 = 11桁
                else:
                    access_nid = f'0{mcode[1:]}'   # "0"   + 11桁コードの下10桁 = 11桁

            # 最近傍道路ノードを検索
            _, idx = tree.query([c_lon, c_lat])
            n_lon = lons[idx]
            n_lat = lats[idx]
            road_nid = node_ids[idx]

            dist_m = max(1, round(haversine_m(c_lon, c_lat, n_lon, n_lat)))
            tim    = max(1, round(dist_m / ACCESS_KMH * 6.0))  # 0.01分単位

            # 2次メッシュ: メッシュコードの先頭6桁
            mesh2 = mcode[:6]

            uid_f = f'A{i+1:015d}'
            uid_r = f'B{i+1:015d}'

            common = f'0,9,0,0,{dist_m},9,2,0,0,0,0,0,00000000000'

            # 正方向（重心→道路ノード）
            csv_out.write(f'{mesh2},{access_nid},{road_nid},{common},1,{uid_f},{tim}\n')
            link_cnt += 1
            # 逆方向（道路ノード→重心）
            csv_out.write(f'{mesh2},{road_nid},{access_nid},{common},2,{uid_r},{tim}\n')
            link_cnt += 1

            # make_net.f用3列CSV
            dist_real = haversine_m(c_lon, c_lat, n_lon, n_lat)
            mn_out.write(f'{access_nid},{road_nid},{dist_real:.6f}\n')

            # GeoJSON
            if not first_feature:
                gjson_out.write(',\n')
            first_feature = False
            gjson_out.write(
                '{"type":"Feature","geometry":{"type":"LineString","coordinates":['
                f'[{c_lon},{c_lat}],[{n_lon},{n_lat}]'
                ']},"properties":{'
                f'"access_node":"{access_nid}",'
                f'"road_node":"{road_nid}",'
                f'"mesh_code":"{mcode}",'
                f'"mode":"{("census" if CENSUS_MODE else "L"+str(LEVEL))}",'
                f'"dist_m":{dist_m},'
                f'"time_001min":{tim}'
                '}}'
            )

            if (i + 1) % 5000 == 0:
                pct = (i + 1) * 100 // total
                print(f'      {i+1:,} / {total:,}  ({pct}%)  リンク数: {link_cnt:,}')

        gjson_out.write('\n]}\n')

    print(f'      アクセスリンク出力数: {link_cnt:,}  ({time.time()-t0:.1f}s)')

    # ── [4/4] GeoParquet出力 ────────────────────────────────────
    import geopandas as gpd
    import os
    OUT_ACCESS_PARQUET = OUT_DIR / f'KSJ_N13-24_{_prefix}_アクセスリンク_{_suffix}.parquet'
    print(f'\n[4/4] GeoParquet出力中 ...')
    gdf = gpd.read_file(OUT_ACCESS_GEOJSON)
    gdf.to_parquet(OUT_ACCESS_PARQUET)
    print(f'      Saved: {OUT_ACCESS_PARQUET.name}')

    # ── 統計 ──────────────────────────────────────────────────
    elapsed = time.time() - t0
    print(f'\n===== 完了 ({elapsed:.0f}s) =====')
    if CENSUS_MODE:
        print(f'  censusメッシュ数（100m）  : {total:,}')
    else:
        mesh_size_str = {3: '1km', 4: '500m', 5: '250m', 6: '125m'}[LEVEL]
        range_label = '全国' if NATIONWIDE_MODE else (str(MESH_LIST) if MULTI_MESH_MODE else MESH_CODE)
        print(f'  対象範囲                  : {range_label}')
        print(f'  {LEVEL}次メッシュ数（{mesh_size_str}）      : {total:,}')
    print(f'  出力アクセスリンク数      : {link_cnt:,}（双方向）')
    print()
    for p in [OUT_ACCESS_CSV, OUT_ACCESS_MAKENET, OUT_ACCESS_GEOJSON, OUT_ACCESS_PARQUET]:
        size_mb = os.path.getsize(p) / 1024 / 1024
        print(f'  {p.name:<55s}  {size_mb:5.1f} MB')


if __name__ == '__main__':
    main()
