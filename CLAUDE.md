# 時刻スライス型到達圏分析 — Claude Code 向けプロジェクト仕様

## プロジェクト概要

「AからBまでT_max分で移動した」実績から、移動中に立ち寄れた可能性のあるエリアを 125m メッシュ（L6）で出力するツール。捜査支援（防犯カメラ確認・聞き込みエリアの絞り込み）を想定。

## ディレクトリ構成

```
src/
  timeslice_search.py   メイン分析スクリプト（双方向Dijkstra）
  ksj_to_network_csv.py 国土数値情報 GeoJSON → 道路リンク/ノード parquet 変換
  make_access_links.py  L6 アクセスリンク（メッシュ重心→最寄道路ノード）生成

data/
  saitama/              埼玉県サンプルネットワーク（同梱済み）
  prefecture.parquet    都道府県境界（--pref クリップ用）
  city.parquet          市区町村境界（--city クリップ用）

input/                  国土数値情報 GeoJSON 配置場所（gitignored）
network/                生成ネットワークデータ（gitignored）
output/                 分析出力（gitignored）
```

## 主要スクリプトの仕様

### timeslice_search.py

- **入力**: 道路リンク / 道路ノード / L6 アクセスリンク（各 parquet）
- **アルゴリズム**: 前向き Dijkstra（A→全ノード）+ 後向き Dijkstra（B→全ノード・G.T）
- **出力**: `output/arrival_map_*.parquet/.qml`、`output/timeslice_*.parquet/.qml`
- **パス定数**: `REPO_ROOT = Path(__file__).parent.parent`（`src/` の 1 つ上がリポジトリルート）
- **デフォルトデータ**: `data/saitama/` の埼玉県サンプル

### ksj_to_network_csv.py

- **入力**: `input/N13-24_{mesh}/N13-24_{mesh}.geojson`
- **出力**: `network/{case}/KSJ_N13-24_{case}_道路リンク.parquet` など
- **パス定数**: `BASE_DIR = Path(__file__).parent.parent`（同上）
- **クリップ**: `data/prefecture.parquet`（--pref）、`data/city.parquet`（--city）
- **重要**: 全道路を双方向リンクとして生成（国土数値情報に一方通行フィールドなし）

### make_access_links.py

- **入力**: `network/{case}/KSJ_N13-24_{case}_道路ノード.csv`（ksj_to_network_csv.py の出力）
- **出力**: `network/{case}/KSJ_N13-24_{case}_アクセスリンク_L6.parquet` など
- **パス定数**: `BASE_DIR = Path(__file__).parent.parent`（同上）

## 新しいエリアで分析する手順

1. 国土数値情報（https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N13-v2_1.html）から対象都道府県の N13-24 GeoJSON をダウンロードし `input/` に配置
2. `python3 src/ksj_to_network_csv.py --meshes {meshes} --case {name} --pref {pref}` でネットワーク生成
3. `python3 src/make_access_links.py --meshes {meshes} --case {name} --level 6 --pref {pref}` でアクセスリンク生成
4. `python3 src/timeslice_search.py --links network/{name}/... --nodes ... --access ... --orig-lat ... --dest-lat ... --tmax {T}` で分析

## 出力ファイルの構造

**timeslice_*.parquet** カラム:
- `mesh_code`: L6 メッシュコード（11 桁）
- `dist_a_min` / `dist_b_min`: 始点 A・終点 B からの最短到達時間（分）
- `t_arrive_earliest` / `t_arrive_latest`: 最早・最遅到着可能時刻（分）
- `slack_min` / `slack_rank`: 余裕時間・ランク（通過不可は NaN / 空文字）
- `in_t10` 〜 `in_tXX`: 各時刻スライスで存在可能か（True/False）
- `geometry`: メッシュポリゴン（EPSG:4326）

## 制約

- **一方通行未考慮**: 国土数値情報に一方通行フィールドが存在しないため全道路双方向
- **メッシュレベル**: L6（125m）固定。変更する場合は `--level` オプションと `timeslice_search.py` のアクセスリンクパスを合わせて変更する
- **速度モデル**: vehicle=道路種別・幅員による速度テーブル、walk=3.6 km/h 一律
