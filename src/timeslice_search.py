#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
時刻スライス型到達圏分析スクリプト

始点A・終点B・移動時間実績 T_max をもとに、各時刻における
存在可能エリアを L6 メッシュ（125m）GeoParquet で出力する。

【アルゴリズム】
  前向きDijkstra : A → 全ノード → dist_A[m]
  後向きDijkstra : B → 全ノード（逆向きグラフ）→ dist_B[m]

  各メッシュ m の判定:
    slack[m] = T_max - dist_A[m] - dist_B[m]
    slack >= 0 → 通過可能
    時刻 t に存在可能: dist_A[m] <= t AND dist_B[m] <= T_max - t

【使い方】
  python3 timeslice_search.py                           # デフォルト（埼玉県庁→東松山市役所 60分）
  python3 timeslice_search.py --tmax 70                 # T_max を変更
  python3 timeslice_search.py --tmax 60,70,80           # 複数 T_max を一括出力（Dijkstra は1回）
  python3 timeslice_search.py --slice 5                 # スライス間隔を5分に
"""

import argparse
import json
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra as sp_dijkstra
from scipy.spatial import KDTree
from shapely.geometry import box

SCRIPT_DIR     = Path(__file__).parent
REPO_ROOT      = SCRIPT_DIR.parent
SAMPLE_DIR     = REPO_ROOT / "network" / "saitama"
DEFAULT_LINKS  = str(SAMPLE_DIR / "KSJ_N13-24_saitama_all_道路リンク.parquet")
DEFAULT_NODES  = str(SAMPLE_DIR / "KSJ_N13-24_saitama_all_道路ノード.parquet")
DEFAULT_ACCESS = str(SAMPLE_DIR / "KSJ_N13-24_saitama_all_アクセスリンク_L6.parquet")
DEFAULT_TMAX   = 60.0
DEFAULT_SLICE  = 10
DEFAULT_OUT    = str(REPO_ROOT / "output")

# 到達時間ランク（10 分刻み・赤→橙→黄→緑→シアン→青→濃紫）
VEHICLE_BINS = list(range(0, 100, 10)) + [float("inf")]  # 0,10,20,...,90,inf
VEHICLE_LABELS = list(range(len(VEHICLE_BINS) - 1))      # 0〜9 の 10 段階
VEHICLE_COLORS = [
    (255,   0,   0),  #  0〜 10分
    (255,  64,   0),  # 10〜 20分
    (255, 128,   0),  # 20〜 30分
    (255, 192,   0),  # 30〜 40分
    (255, 255,   0),  # 40〜 50分
    (192, 255,   0),  # 50〜 60分
    (  0, 204,   0),  # 60〜 70分
    (  0, 204, 128),  # 70〜 80分
    (  0, 204, 204),  # 80〜 90分
    ( 68,   0,  85),  # 90分超
]

# slack 10 分刻み色テーブル（赤→橙→黄→黄緑→緑→青）
SLACK_COLORS = [
    (215,  25,  28),  # 0〜10分  赤（ほぼ必須経路）
    (253, 174,  97),  # 10〜20分 橙
    (255, 255, 191),  # 20〜30分 黄
    (166, 217, 106),  # 30〜40分 黄緑
    ( 26, 152,  80),  # 40〜50分 緑
    (  0,  85, 255),  # 50分超   青
]


# ──────────────────────────────────────────────
# L6 メッシュポリゴン生成
# ──────────────────────────────────────────────
def compute_l6_polygons(codes_series: pd.Series) -> gpd.GeoDataFrame:
    codes = codes_series.astype(str).str.zfill(11)
    pp = codes.str[0:2].astype(int).values
    qq = codes.str[2:4].astype(int).values
    r  = codes.str[4].astype(int).values
    s  = codes.str[5].astype(int).values
    t  = codes.str[6].astype(int).values
    u  = codes.str[7].astype(int).values
    v  = codes.str[8].astype(int).values
    w  = codes.str[9].astype(int).values
    x  = codes.str[10].astype(int).values

    lat = pp / 1.5
    lon = (qq + 100).astype(float)
    dlat2 = (2.0/3.0)/8;  dlon2 = 1.0/8
    lat += r * dlat2;     lon += s * dlon2
    dlat3 = dlat2/10;     dlon3 = dlon2/10
    lat += t * dlat3;     lon += u * dlon3
    dlat4 = dlat3/2;      dlon4 = dlon3/2
    lat += ((v-1)//2) * dlat4;  lon += ((v-1)%2) * dlon4
    dlat5 = dlat4/2;      dlon5 = dlon4/2
    lat += ((w-1)//2) * dlat5;  lon += ((w-1)%2) * dlon5
    dlat6 = dlat5/2;      dlon6 = dlon5/2
    lat += ((x-1)//2) * dlat6;  lon += ((x-1)%2) * dlon6

    geoms = [box(lo, la, lo+dlon6, la+dlat6) for lo, la in zip(lon, lat)]
    return gpd.GeoDataFrame({"mesh_code": codes.values}, geometry=geoms, crs="EPSG:4326")


# ──────────────────────────────────────────────
# グラフ構築（全道路双方向・KSJに一方通行フィールドなし）
# ──────────────────────────────────────────────
def build_sparse_graph(links: gpd.GeoDataFrame):
    n1 = links["node1"].astype(int).to_numpy()
    n2 = links["node2"].astype(int).to_numpy()
    ws = (links["time_001min"].astype(float) * 0.01).to_numpy(dtype=np.float64)

    src = np.concatenate([n1, n2])
    dst = np.concatenate([n2, n1])
    w   = np.concatenate([ws, ws])

    unique = np.unique(np.concatenate([src, dst]))
    n2i    = {int(n): i for i, n in enumerate(unique.tolist())}
    n_v    = len(unique)

    rows = np.array([n2i[int(n)] for n in src], dtype=np.int32)
    cols = np.array([n2i[int(n)] for n in dst], dtype=np.int32)
    G    = csr_matrix((w, (rows, cols)), shape=(n_v, n_v))
    return unique, n2i, G


# ──────────────────────────────────────────────
# 最近傍道路ノード検索
# ──────────────────────────────────────────────
def nearest_road_node(nodes: gpd.GeoDataFrame, lat: float, lon: float):
    coords = np.array([[p.y, p.x] for p in nodes.geometry])
    dist, i = KDTree(coords).query([lat, lon])
    row = nodes.iloc[i]
    return int(row["node_id"]), float(row.geometry.y), float(row.geometry.x), float(dist * 111000)


# ──────────────────────────────────────────────
# 到達時間マップ出力（A・B 単独）
# ──────────────────────────────────────────────
def write_arrival_qml(qml_path: Path, label_prefix: str) -> None:
    n = len(VEHICLE_LABELS)

    def label_text(i):
        lo = int(VEHICLE_BINS[i])
        hi = VEHICLE_BINS[i + 1]
        return f"{lo}分超" if hi == float("inf") else f"{lo}〜{int(hi)}分"

    cats = "\n".join(
        f'      <category symbol="{i}" value="{VEHICLE_LABELS[i]}" '
        f'label="{label_text(i)}" render="true"/>'
        for i in range(n)
    )
    cats += f'\n      <category symbol="{n}" value="" label="到達不能" render="true"/>'

    syms = []
    for i, (r, g, b) in enumerate(VEHICLE_COLORS):
        syms.append(
            f'      <symbol name="{i}" type="fill" alpha="0.75" '
            f'clip_to_extent="1" is_animated="0" frame_rate="10">\n'
            f'        <data_defined_properties><Option type="Map">'
            f'<Option name="name" type="QString" value=""/>'
            f'<Option name="properties"/>'
            f'<Option name="type" type="QString" value="collection"/>'
            f'</Option></data_defined_properties>\n'
            f'        <layer class="SimpleFill" enabled="1" pass="0" locked="0">\n'
            f'          <Option type="Map">\n'
            f'            <Option name="color" type="QString" value="{r},{g},{b},255"/>\n'
            f'            <Option name="outline_style" type="QString" value="no"/>\n'
            f'            <Option name="style" type="QString" value="solid"/>\n'
            f'          </Option>\n'
            f'        </layer>\n'
            f'      </symbol>'
        )
    syms.append(
        f'      <symbol name="{n}" type="fill" alpha="0.75" '
        f'clip_to_extent="1" is_animated="0" frame_rate="10">\n'
        f'        <data_defined_properties><Option type="Map">'
        f'<Option name="name" type="QString" value=""/>'
        f'<Option name="properties"/>'
        f'<Option name="type" type="QString" value="collection"/>'
        f'</Option></data_defined_properties>\n'
        f'        <layer class="SimpleFill" enabled="1" pass="0" locked="0">\n'
        f'          <Option type="Map">\n'
        f'            <Option name="color" type="QString" value="170,170,170,180"/>\n'
        f'            <Option name="outline_style" type="QString" value="no"/>\n'
        f'            <Option name="style" type="QString" value="solid"/>\n'
        f'          </Option>\n'
        f'        </layer>\n'
        f'      </symbol>'
    )

    content = (
        "<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>\n"
        '<qgis version="3.34.0" styleCategories="Symbology">\n'
        '  <renderer-v2 type="categorizedSymbol" attr="到達時間ランク" '
        'forceraster="0" symbollevels="0" usingSymbolLevels="0" enableorderby="0">\n'
        '    <categories>\n'
        f'{cats}\n'
        '    </categories>\n'
        '    <symbols>\n'
        + "\n".join(syms) + "\n"
        '    </symbols>\n'
        '    <rotation/>\n'
        '    <sizescale/>\n'
        '  </renderer-v2>\n'
        '  <blendMode>0</blendMode>\n'
        '  <featureBlendMode>0</featureBlendMode>\n'
        '  <layerOpacity>1</layerOpacity>\n'
        '</qgis>\n'
    )
    qml_path.write_text(content, encoding="utf-8")


def save_arrival_map(mesh_codes, dist_mesh, out_parquet: Path, out_qml: Path,
                     label: str) -> None:
    reachable = np.isfinite(dist_mesh)
    dist_safe = np.where(reachable, dist_mesh, 0.0)
    rank_idx  = np.minimum(np.digitize(dist_safe, VEHICLE_BINS[1:]), len(VEHICLE_LABELS) - 1)
    rank_str  = np.where(reachable, rank_idx.astype(str), "")

    df = pd.DataFrame({
        "mesh_code":    mesh_codes,
        "到達時間_min": np.where(reachable, np.round(dist_mesh, 2), np.nan),
        "到達時間ランク": rank_str,
    })

    mesh_gdf = compute_l6_polygons(pd.Series(mesh_codes))
    gdf = mesh_gdf.merge(df, on="mesh_code").set_crs("EPSG:4326")
    gdf.to_parquet(out_parquet)
    write_arrival_qml(out_qml, label)

    reachable_cnt = int(reachable.sum())
    sz = out_parquet.stat().st_size // 1024
    print(f"→ {out_parquet.name}  ({reachable_cnt:,}メッシュ到達可能, {sz}KB)")
    print(f"→ {out_qml.name}")


# ──────────────────────────────────────────────
# QML 生成（通過可能=赤・通過不可=グレーの2色）
# ──────────────────────────────────────────────
def write_qml(qml_path: Path, tmax: float, interval: int) -> None:
    n_col = len(SLACK_COLORS)
    passable_vals = [str(i * interval) for i in range(n_col)]

    cats = "\n".join(
        f'      <category symbol="0" value="{v}" label="通過可能" render="true"/>'
        for v in passable_vals
    )
    cats += '\n      <category symbol="1" value="" label="通過不可" render="true"/>'

    def sym(name, r, g, b, a):
        return (
            f'      <symbol name="{name}" type="fill" alpha="0.75" '
            f'clip_to_extent="1" is_animated="0" frame_rate="10">\n'
            f'        <data_defined_properties><Option type="Map">'
            f'<Option name="name" type="QString" value=""/>'
            f'<Option name="properties"/>'
            f'<Option name="type" type="QString" value="collection"/>'
            f'</Option></data_defined_properties>\n'
            f'        <layer class="SimpleFill" enabled="1" pass="0" locked="0">\n'
            f'          <Option type="Map">\n'
            f'            <Option name="color" type="QString" value="{r},{g},{b},{a}"/>\n'
            f'            <Option name="outline_style" type="QString" value="no"/>\n'
            f'            <Option name="style" type="QString" value="solid"/>\n'
            f'          </Option>\n'
            f'        </layer>\n'
            f'      </symbol>'
        )

    syms = [
        sym("0", 255, 0, 0, 255),        # 通過可能: 赤
        sym("1", 170, 170, 170, 180),     # 通過不可: グレー
    ]

    content = (
        "<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>\n"
        '<qgis version="3.34.0" styleCategories="Symbology">\n'
        '  <renderer-v2 type="categorizedSymbol" attr="slack_rank" '
        'forceraster="0" symbollevels="0" usingSymbolLevels="0" enableorderby="0">\n'
        '    <categories>\n'
        f'{cats}\n'
        '    </categories>\n'
        '    <symbols>\n'
        + "\n".join(syms) + "\n"
        '    </symbols>\n'
        '    <rotation/>\n'
        '    <sizescale/>\n'
        '  </renderer-v2>\n'
        '  <blendMode>0</blendMode>\n'
        '  <featureBlendMode>0</featureBlendMode>\n'
        '  <layerOpacity>1</layerOpacity>\n'
        '</qgis>\n'
    )
    qml_path.write_text(content, encoding="utf-8")


def write_od_qml(qml_path: Path) -> None:
    content = (
        '<!DOCTYPE qgis PUBLIC \'http://mrcc.com/qgis.dtd\' \'SYSTEM\'>\n'
        '<qgis version="3.34.0" styleCategories="Symbology|Labeling">\n'
        '  <renderer-v2 type="singleSymbol" forceraster="0" symbollevels="0"'
        ' usingSymbolLevels="0" enableorderby="0">\n'
        '    <symbols>\n'
        '      <symbol name="0" type="marker" alpha="1" clip_to_extent="1"'
        ' is_animated="0" frame_rate="10">\n'
        '        <data_defined_properties><Option type="Map">'
        '<Option name="name" type="QString" value=""/>'
        '<Option name="properties"/>'
        '<Option name="type" type="QString" value="collection"/>'
        '</Option></data_defined_properties>\n'
        '        <layer class="SimpleMarker" enabled="1" pass="0" locked="0">\n'
        '          <Option type="Map">\n'
        '            <Option name="color" type="QString" value="255,255,255,255"/>\n'
        '            <Option name="outline_color" type="QString" value="0,0,0,255"/>\n'
        '            <Option name="outline_width" type="QString" value="0.4"/>\n'
        '            <Option name="outline_width_unit" type="QString" value="MM"/>\n'
        '            <Option name="size" type="QString" value="3"/>\n'
        '            <Option name="size_unit" type="QString" value="MM"/>\n'
        '            <Option name="name" type="QString" value="circle"/>\n'
        '          </Option>\n'
        '        </layer>\n'
        '      </symbol>\n'
        '    </symbols>\n'
        '    <rotation/>\n'
        '    <sizescale/>\n'
        '  </renderer-v2>\n'
        '  <labeling type="simple">\n'
        '    <settings calloutType="simple">\n'
        '      <text-style fieldName="name" fontFamily="sans-serif" fontSize="10"'
        ' fontWeight="75" textColor="0,0,0,255" namedStyle="Bold"'
        ' textOpacity="1" blendMode="0" isExpression="0"'
        ' fontLetterSpacing="0" fontWordSpacing="0"'
        ' fontUnderline="0" fontStrikeout="0" fontItalic="0"'
        ' fontSizeUnit="Point" fontSizeMapUnitScale="3x:0,0,0,0,0,0">\n'
        '        <text-buffer bufferDraw="1" bufferSize="1" bufferSizeUnits="MM"'
        ' bufferColor="255,255,255,255" bufferOpacity="1" bufferBlendMode="0"'
        ' bufferNoFill="0" bufferJoinStyle="128"/>\n'
        '        <background shapeDraw="0"/>\n'
        '        <shadow shadowDraw="0"/>\n'
        '      </text-style>\n'
        '      <text-format autoWrapLength="0" useMaxLineLengthForAutoWrap="1"'
        ' addDirectionSymbol="0" leftDirectionSymbol="&lt;"'
        ' rightDirectionSymbol="&gt;" reverseDirectionSymbol="0"'
        ' placeDirectionSymbol="0" formatNumbers="0" decimals="3"'
        ' plusSign="0" multilineAlign="3"/>\n'
        '      <placement placement="2" offsetType="0" xOffset="0" yOffset="0"'
        ' offsetUnits="MM" dist="2" distUnits="MM" distMapUnitScale="3x:0,0,0,0,0,0"'
        ' repeatDistance="0" repeatDistanceUnits="MM"'
        ' repeatDistanceMapUnitScale="3x:0,0,0,0,0,0"'
        ' maxCurvedCharAngleIn="25" maxCurvedCharAngleOut="-25"'
        ' priority="5" predefinedPositionOrder="TR,TL,BR,BL,R,L,TSR,BSR"'
        ' fitInPolygonOnly="0" overrunDistance="0" overrunDistanceUnit="MM"'
        ' overrunDistanceMapUnitScale="3x:0,0,0,0,0,0"'
        ' labelOffsetMapUnitScale="3x:0,0,0,0,0,0"'
        ' polygonPlacementFlags="2" allowDegraded="0"'
        ' geometryGenerator="" geometryGeneratorEnabled="0"'
        ' geometryGeneratorType="PointGeometry" layerType="PointGeometry"'
        ' centroidWhole="0" centroidInside="0"'
        ' overlapHandling="PreventOverlap" zIndex="0"/>\n'
        '      <rendering obstacle="1" obstacleFactor="1" obstacleType="1"'
        ' scaleVisibility="0" minScale="1" maxScale="0"'
        ' limitNumLabels="0" maxNumLabels="2000"'
        ' displayAll="0" upsidedownLabels="0"'
        ' fontMinPixelSize="3" fontMaxPixelSize="10000"'
        ' mergeLines="0" drawLabels="1" labelPerPart="0"'
        ' scaleMin="0" scaleMax="0"/>\n'
        '      <dd_properties>\n'
        '        <Option type="Map"><Option name="name" type="QString" value=""/>'
        '<Option name="properties"/>'
        '<Option name="type" type="QString" value="collection"/></Option>\n'
        '      </dd_properties>\n'
        '      <callout type="simple">\n'
        '        <Option type="Map"><Option name="anchorPoint" type="QString" value="pole_of_inaccessibility"/>'
        '<Option name="blendMode" type="int" value="0"/>'
        '<Option name="ddProperties" type="Map">'
        '<Option name="name" type="QString" value=""/>'
        '<Option name="properties"/>'
        '<Option name="type" type="QString" value="collection"/>'
        '</Option>'
        '<Option name="drawToAllParts" type="bool" value="false"/>'
        '<Option name="enabled" type="QString" value="0"/>'
        '<Option name="labelAnchorPoint" type="QString" value="point_on_exterior"/>'
        '<Option name="lineSymbol" type="QString" value="&lt;symbol name=&quot;_&quot; type=&quot;line&quot; alpha=&quot;1&quot; clip_to_extent=&quot;1&quot; is_animated=&quot;0&quot; frame_rate=&quot;10&quot;&gt;&lt;data_defined_properties&gt;&lt;Option type=&quot;Map&quot;&gt;&lt;Option name=&quot;name&quot; type=&quot;QString&quot; value=&quot;&quot;/&gt;&lt;Option name=&quot;properties&quot;/&gt;&lt;Option name=&quot;type&quot; type=&quot;QString&quot; value=&quot;collection&quot;/&gt;&lt;/Option&gt;&lt;/data_defined_properties&gt;&lt;layer class=&quot;SimpleLine&quot; enabled=&quot;1&quot; pass=&quot;0&quot; locked=&quot;0&quot;&gt;&lt;Option type=&quot;Map&quot;&gt;&lt;Option name=&quot;line_color&quot; type=&quot;QString&quot; value=&quot;60,60,60,255&quot;/&gt;&lt;Option name=&quot;line_width&quot; type=&quot;QString&quot; value=&quot;0.3&quot;/&gt;&lt;/Option&gt;&lt;/layer&gt;&lt;/symbol&gt;"/>'
        '<Option name="minLength" type="double" value="0"/>'
        '<Option name="minLengthMapUnitScale" type="QString" value="3x:0,0,0,0,0,0"/>'
        '<Option name="minLengthUnit" type="QString" value="MM"/>'
        '<Option name="offsetFromAnchor" type="double" value="0"/>'
        '<Option name="offsetFromAnchorMapUnitScale" type="QString" value="3x:0,0,0,0,0,0"/>'
        '<Option name="offsetFromAnchorUnit" type="QString" value="MM"/>'
        '<Option name="offsetFromLabel" type="double" value="0"/>'
        '<Option name="offsetFromLabelMapUnitScale" type="QString" value="3x:0,0,0,0,0,0"/>'
        '<Option name="offsetFromLabelUnit" type="QString" value="MM"/>'
        '</Option>\n'
        '      </callout>\n'
        '    </settings>\n'
        '  </labeling>\n'
        '  <blendMode>0</blendMode>\n'
        '  <featureBlendMode>0</featureBlendMode>\n'
        '  <layerOpacity>1</layerOpacity>\n'
        '</qgis>\n'
    )
    qml_path.write_text(content, encoding="utf-8")


# ──────────────────────────────────────────────
# メイン
# ──────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description="時刻スライス型到達圏分析")
    ap.add_argument("--links",    default=DEFAULT_LINKS,  help="道路リンク parquet")
    ap.add_argument("--nodes",    default=DEFAULT_NODES,  help="道路ノード parquet")
    ap.add_argument("--access",   default=DEFAULT_ACCESS, help="L6アクセスリンク parquet")
    ap.add_argument("--orig-lat", type=float, default=35.8578, help="始点緯度（埼玉県庁）")
    ap.add_argument("--orig-lon", type=float, default=139.6490, help="始点経度（埼玉県庁）")
    ap.add_argument("--dest-lat", type=float, default=36.0420, help="終点緯度（東松山市役所）")
    ap.add_argument("--dest-lon", type=float, default=139.4006, help="終点経度（東松山市役所）")
    ap.add_argument("--tmax",      type=str,   default=str(DEFAULT_TMAX),
                    help="移動時間実績 T_max（分）。カンマ区切りで複数指定可（例: 60,70,80）")
    ap.add_argument("--slice",     type=int,   default=DEFAULT_SLICE, help="時刻スライス間隔（分）")
    ap.add_argument("--orig-name", default="埼玉県庁",   help="始点名（出力ファイル名に使用）")
    ap.add_argument("--dest-name", default="東松山市役所", help="終点名（出力ファイル名に使用）")
    ap.add_argument("--out-dir",   default=DEFAULT_OUT, help="出力ディレクトリ")
    args = ap.parse_args()
    tmax_raw = [float(t) for t in args.tmax.replace(" ", "").split(",") if t]

    t0      = time.time()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── データ読み込み ────────────────────────
    print("道路リンク読み込み中...")
    links = gpd.read_parquet(args.links)
    print(f"  {len(links):,} 本")

    print("道路ノード読み込み中...")
    nodes = gpd.read_parquet(args.nodes)
    print(f"  {len(nodes):,} 件")

    print("アクセスリンク読み込み中...")
    access = gpd.read_parquet(args.access)
    print(f"  {len(access):,} 件  ({time.time()-t0:.1f}s)")

    # ── グラフ構築 ────────────────────────────
    print("グラフ構築中...")
    _, n2i, G = build_sparse_graph(links)
    print(f"  ノード: {G.shape[0]:,}  エッジ: {G.nnz:,}  ({time.time()-t0:.1f}s)")

    # ── OD ノード検索 ─────────────────────────
    orig_nid, o_lat, o_lon, o_snap = nearest_road_node(nodes, args.orig_lat, args.orig_lon)
    dest_nid, d_lat, d_lon, d_snap = nearest_road_node(nodes, args.dest_lat, args.dest_lon)
    print(f"始点ノード: {orig_nid}  ({o_lat:.5f},{o_lon:.5f})  snap={o_snap:.0f}m")
    print(f"終点ノード: {dest_nid}  ({d_lat:.5f},{d_lon:.5f})  snap={d_snap:.0f}m")

    if orig_nid not in n2i:
        raise ValueError(f"始点ノード {orig_nid} がグラフに存在しません")
    if dest_nid not in n2i:
        raise ValueError(f"終点ノード {dest_nid} がグラフに存在しません")

    orig_idx = n2i[orig_nid]
    dest_idx = n2i[dest_nid]

    # ── 双方向ダイクストラ ────────────────────
    print("前向きDijkstra（A → 全ノード）...")
    dist_A = sp_dijkstra(G,   directed=True, indices=orig_idx)
    print(f"  完了  ({time.time()-t0:.1f}s)")

    print("後向きDijkstra（B → 全ノード・逆向きグラフ）...")
    dist_B = sp_dijkstra(G.T, directed=True, indices=dest_idx)
    print(f"  完了  ({time.time()-t0:.1f}s)")

    t_min_road = dist_A[dest_idx]
    print(f"最短経路時間: {t_min_road:.1f}分")
    tmax_list = sorted({t for t in tmax_raw if t >= t_min_road})
    if not tmax_list:
        raise ValueError(f"すべての T_max が最短経路時間({t_min_road:.1f}分)より小さいです")
    print(f"出力対象 T_max: {tmax_list}")

    # ── メッシュ到達時間計算（ベクトル化・tmax 非依存）────
    print("メッシュ別到達時間計算中...")
    road_nids  = access["road_node"].astype(int).to_numpy()
    acc_time   = (access["time_001min"].astype(float) * 0.01).to_numpy()
    mesh_codes = access["mesh_code"].astype(str).to_numpy()

    graph_idxs = np.array([n2i.get(int(rn), -1) for rn in road_nids], dtype=np.int32)
    valid      = graph_idxs >= 0
    safe_idxs  = np.where(valid, graph_idxs, 0)
    da_road    = np.where(valid, dist_A[safe_idxs], np.inf)
    db_road    = np.where(valid, dist_B[safe_idxs], np.inf)

    dist_A_mesh = np.where(da_road < np.inf, da_road + acc_time, np.inf)
    dist_B_mesh = np.where(db_road < np.inf, db_road + acc_time, np.inf)
    print(f"  ({time.time()-t0:.1f}s)")

    # ── 到達時間マップ（A・B 単独）出力 ──────
    print("\n到達時間マップ（単独）出力中...")
    save_arrival_map(
        mesh_codes, dist_A_mesh,
        out_dir / f"arrival_map_{args.orig_name}.parquet",
        out_dir / f"arrival_map_{args.orig_name}.qml",
        args.orig_name,
    )
    save_arrival_map(
        mesh_codes, dist_B_mesh,
        out_dir / f"arrival_map_{args.dest_name}.parquet",
        out_dir / f"arrival_map_{args.dest_name}.qml",
        args.dest_name,
    )
    print(f"  ({time.time()-t0:.1f}s)\n")

    # ── メッシュポリゴン（一度だけ生成）────────
    print("メッシュポリゴン生成中...")
    mesh_gdf = compute_l6_polygons(pd.Series(mesh_codes))
    print(f"  ({time.time()-t0:.1f}s)\n")

    # ── tmax ごとに時刻スライス出力 ────────────
    n_col    = len(SLACK_COLORS)
    bin_strs = np.array([str(i * args.slice) for i in range(n_col)])

    for tmax in tmax_list:
        print(f"── T_max={int(tmax)}分 ──────────────────────────────")
        slack    = tmax - dist_A_mesh - dist_B_mesh
        feasible = slack >= 0
        print(f"  通過可能: {feasible.sum():,} / {len(access):,} メッシュ  余裕: {tmax - t_min_road:.1f}分")

        slices = list(range(args.slice, int(tmax), args.slice))
        slice_data = {}
        for t in slices:
            col  = f"in_t{t:02d}"
            in_s = feasible & (dist_A_mesh <= t) & (dist_B_mesh <= tmax - t)
            slice_data[col] = in_s
            print(f"  t={t:2d}分: {in_s.sum():,} メッシュ")

        slack_safe   = np.where(feasible, slack, 0.0)
        s_rank_idx   = np.where(feasible,
                                np.minimum((slack_safe / args.slice).astype(int), n_col - 1),
                                -1)
        slack_rank   = np.where(s_rank_idx >= 0, bin_strs[np.clip(s_rank_idx, 0, n_col - 1)], "")

        da_safe      = np.where(feasible, dist_A_mesh, 0.0)
        da_rank_idx  = np.minimum(np.digitize(da_safe, VEHICLE_BINS[1:]), len(VEHICLE_LABELS) - 1)
        dist_a_rank  = np.where(feasible, da_rank_idx.astype(str), "")

        df = pd.DataFrame({
            "mesh_code":         mesh_codes,
            "dist_a_min":        np.where(feasible, np.round(dist_A_mesh, 2), np.nan),
            "dist_a_rank":       dist_a_rank,
            "dist_b_min":        np.where(feasible, np.round(dist_B_mesh, 2), np.nan),
            "t_arrive_earliest": np.where(feasible, np.round(dist_A_mesh, 2), np.nan),
            "t_arrive_latest":   np.where(feasible, np.round(tmax - dist_B_mesh, 2), np.nan),
            "slack_min":         np.where(feasible, np.round(slack, 2), np.nan),
            "slack_rank":        slack_rank,
            **{col: v for col, v in slice_data.items()},
        })

        gdf = mesh_gdf.merge(df, on="mesh_code").set_crs("EPSG:4326")

        stem        = f"timeslice_{args.orig_name}_{args.dest_name}_tmax{int(tmax)}"
        out_parquet = out_dir / f"{stem}.parquet"
        out_qml     = out_dir / f"{stem}.qml"

        gdf.to_parquet(out_parquet)
        write_qml(out_qml, tmax, args.slice)

        sz = out_parquet.stat().st_size // 1024
        print(f"→ {out_parquet.name}  ({feasible.sum():,}メッシュ, {sz}KB)")
        print(f"→ {out_qml.name}\n")

    # ── OD ポイント出力 ──────────────────────────
    od_geojson = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": [args.orig_lon, args.orig_lat]},
                "properties": {
                    "type": "origin",
                    "name": args.orig_name,
                    "lat": args.orig_lat,
                    "lon": args.orig_lon,
                    "snap_lat": o_lat,
                    "snap_lon": o_lon,
                    "snap_m": round(o_snap),
                },
            },
            {
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": [args.dest_lon, args.dest_lat]},
                "properties": {
                    "type": "destination",
                    "name": args.dest_name,
                    "lat": args.dest_lat,
                    "lon": args.dest_lon,
                    "snap_lat": d_lat,
                    "snap_lon": d_lon,
                    "snap_m": round(d_snap),
                },
            },
        ],
    }
    od_path = out_dir / f"od_points_{args.orig_name}_{args.dest_name}.geojson"
    with open(od_path, "w", encoding="utf-8") as f:
        json.dump(od_geojson, f, ensure_ascii=False, indent=2)
    print(f"→ {od_path.name}")

    od_qml_path = out_dir / f"od_points_{args.orig_name}_{args.dest_name}.qml"
    write_od_qml(od_qml_path)
    print(f"→ {od_qml_path.name}")

    print(f"\n総処理時間: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
