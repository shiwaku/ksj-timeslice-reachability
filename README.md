# 時刻スライス型到達圏分析

「AからBまでT_max分で移動した」という実績をもとに、移動途中で立ち寄れた可能性のあるエリアを 125m メッシュ（L6）で出力するツール。防犯カメラ確認・聞き込みエリアの絞り込みなど、捜査支援を想定している。

---

## アルゴリズム

双方向ダイクストラにより各 L6 メッシュの通過可能性を判定する。

```
slack[m] = T_max - dist_A[m] - dist_B[m]

slack[m] >= 0  → 通過可能
slack[m] <  0  → 通過不可

時刻 t に存在可能: dist_A[m] <= t  かつ  dist_B[m] <= T_max - t
```

逆向きグラフを `G.T`（scipy 転置行列）で表現し、`scipy.sparse.csgraph.dijkstra` で実装。外部プログラム不使用。

---

## 制約・注意事項

> ⚠️ **一方通行は考慮されていない**  
> 国土数値情報の道路データには一方通行フィールドが存在しない。このため、全道路を双方向リンクとして扱っており、実際の交通規制は反映されていない。より厳密な分析には OpenStreetMap（`oneway` タグあり）等の別データソースへの切り替えが必要。

---

## ディレクトリ構成

```
（リポジトリルート）
├── README.md
├── CLAUDE.md                     # Claude Code 向けプロジェクト仕様
├── .claude/
│   └── commands/                 # カスタムスラッシュコマンド
│       ├── generate-network.md   # /generate-network: ネットワーク生成
│       └── analyze.md            # /analyze: 時刻スライス分析実行
├── src/
│   ├── timeslice_search.py       # 時刻スライス分析（メイン）
│   ├── ksj_to_network_csv.py     # 国土数値情報 → 道路リンク・ノード parquet
│   └── make_access_links.py      # L6 アクセスリンク生成
├── data/
│   ├── saitama/                  # サンプルデータ（埼玉県・同梱済み）
│   │   ├── KSJ_N13-24_saitama_all_道路リンク.parquet
│   │   ├── KSJ_N13-24_saitama_all_道路ノード.parquet
│   │   └── KSJ_N13-24_saitama_all_アクセスリンク_L6.parquet
│   ├── prefecture.parquet        # 都道府県境界ポリゴン（--pref クリップ用）
│   └── city.parquet              # 市区町村境界ポリゴン（--city クリップ用）
├── input/                        # 国土数値情報 GeoJSON 配置場所（gitignored・各自用意）
├── network/                      # 生成したネットワークデータ（gitignored・再生成可能）
└── output/                       # 分析出力（gitignored・実行時に自動生成）
```

---

## 必要環境

- Python 3.9 以上
- 依存ライブラリ

```bash
pip install geopandas pyarrow scipy
```

---

## サンプルデータでの実行（埼玉県）

`data/saitama/` に埼玉県全道路ネットワークデータ（サンプル）を同梱している。

| ファイル | 内容 | サイズ |
|---|---|---|
| `KSJ_N13-24_saitama_all_道路リンク.parquet` | 道路リンク（949,637本） | 61 MB |
| `KSJ_N13-24_saitama_all_道路ノード.parquet` | 道路ノード（706,418件） | 18 MB |
| `KSJ_N13-24_saitama_all_アクセスリンク_L6.parquet` | L6 アクセスリンク（233,233件） | 7.7 MB |

出典: 国土数値情報 道路データ（N13-24）/ 国土交通省

### デフォルト実行（埼玉県庁 → 東松山市役所・T_max=60分）

```bash
python3 src/timeslice_search.py
```

### T_max を変更

```bash
# T_max を 65 分に変更
python3 src/timeslice_search.py --tmax 65

# 複数 T_max を一括出力（Dijkstra は1回のみ実行）
python3 src/timeslice_search.py --tmax 60,65,70,80
```

### 任意の始点・終点を指定

```bash
python3 src/timeslice_search.py \
  --orig-lat 35.8578 --orig-lon 139.6490 --orig-name 埼玉県庁 \
  --dest-lat 36.0420 --dest-lon 139.4006 --dest-name 東松山市役所 \
  --tmax 60
```

緯度・経度は Google マップで地点を右クリックするとコピーできる。

### 動作確認済み実績（埼玉県庁 → 東松山市役所）

| T_max | 通過可能メッシュ | 処理時間 |
|---|---|---|
| 60 分 | 5,944 件 / 233,233 件 | 約 25 秒 |
| 65 分 | 20,840 件 / 233,233 件 | 〃 |
| 70 分 | 36,134 件 / 233,233 件 | 〃 |

- 最短経路時間: 53.4 分
- グラフ規模: ノード 707,538 / エッジ 1,887,715

---

## 他の都道府県のネットワークデータ作成

埼玉県以外の地域で分析するには、国土数値情報からネットワークデータを生成する必要がある。

### ステップ 1: 国土数値情報のダウンロード

国土交通省「国土数値情報ダウンロードサービス」から道路データ（N13）をダウンロードする。

- **ダウンロード先**: https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N13-v2_1.html
- **対象年度**: 2024 年度版（N13-24）を推奨
- **形式**: GeoJSON

ダウンロードした ZIP を展開し、GeoJSON ファイルを以下の形式で `input/` に配置する。

```
input/
├── N13-24_5133/
│   └── N13-24_5133.geojson
├── N13-24_5134/
│   └── N13-24_5134.geojson
└── ...
```

> **1次メッシュコードの調べ方**  
> 1次メッシュ（約80km×80km 格子）のコードは地域によって異なる。下記の主要都道府県コード表を参考に、対象エリアをカバーするメッシュを選択する。

#### 主要都道府県の 1 次メッシュコード表

| 都道府県 | 1 次メッシュコード |
|---|---|
| 北海道（札幌周辺） | 6441, 6442, 6443, 6444, 6541, 6542, 6543, 6544 |
| 宮城県 | 5640, 5641, 5740, 5741 |
| 東京都 | 5338, 5339, 5438, 5439 |
| 神奈川県 | 5238, 5239, 5338, 5339 |
| 埼玉県 | 5338, 5339, 5438, 5439 |
| 千葉県 | 5239, 5240, 5339, 5340, 5439, 5440 |
| 新潟県 | 5538, 5539, 5638, 5639, 5738, 5739 |
| 愛知県 | 5236, 5237, 5336, 5337, 5436, 5437 |
| 大阪府 | 5135, 5235 |
| 兵庫県 | 5134, 5135, 5234, 5235 |
| 京都府 | 5235, 5335 |
| 広島県 | 5132, 5133, 5232, 5233 |
| 香川県 | 5133, 5134 |
| 鳥取県 | 5233, 5234, 5333, 5334 |
| 福岡県 | 4930, 5030, 5031, 5032 |

> 上記は目安。正確なメッシュリストは国土地理院「地図・空中写真閲覧サービス」等で確認すること。

### ステップ 2: ネットワークデータ生成

リポジトリルートで以下を実行する。

#### 都道府県単位（例: 埼玉県）

```bash
# 道路リンク・ノード生成（全道路・フィルターなし）
python3 src/ksj_to_network_csv.py \
  --meshes 5338,5339,5438,5439 \
  --case saitama_pref \
  --pref 埼玉県

# L6 アクセスリンク生成
python3 src/make_access_links.py \
  --meshes 5338,5339,5438,5439 \
  --case saitama_pref \
  --level 6 \
  --pref 埼玉県
```

出力先: `network/saitama_pref/`

| ファイル | 内容 |
|---|---|
| `KSJ_N13-24_saitama_pref_道路リンク.parquet` | 道路リンク |
| `KSJ_N13-24_saitama_pref_道路ノード.parquet` | 道路ノード |
| `KSJ_N13-24_saitama_pref_アクセスリンク_L6.parquet` | L6 アクセスリンク |

#### 市区町村単位（例: さいたま市）

```bash
python3 src/ksj_to_network_csv.py \
  --meshes 5338,5339,5438,5439 \
  --case saitama_city \
  --city さいたま市 \
  --pref 埼玉県

python3 src/make_access_links.py \
  --meshes 5338,5339,5438,5439 \
  --case saitama_city \
  --level 6 \
  --city さいたま市 \
  --pref 埼玉県
```

### ステップ 3: 生成したデータで分析

```bash
python3 src/timeslice_search.py \
  --links network/saitama_pref/KSJ_N13-24_saitama_pref_道路リンク.parquet \
  --nodes network/saitama_pref/KSJ_N13-24_saitama_pref_道路ノード.parquet \
  --access network/saitama_pref/KSJ_N13-24_saitama_pref_アクセスリンク_L6.parquet \
  --orig-lat 35.8578 --orig-lon 139.6490 --orig-name 埼玉県庁 \
  --dest-lat 36.0420 --dest-lon 139.4006 --dest-name 東松山市役所 \
  --tmax 60
```

### `ksj_to_network_csv.py` 主なオプション

| オプション | 説明 |
|---|---|
| `--meshes 5338,5339,...` | 対象 1 次メッシュコード（カンマ区切り） |
| `--case {name}` | 出力ケース名（`network/{name}/` に出力される） |
| `--pref {都道府県名}` | 都道府県でクリップ（例: `埼玉県`） |
| `--city {市区町村名}` | 市区町村でクリップ（例: `横浜市`） |
| `--mode walk` | 徒歩モード（速度 3.6 km/h。省略時は vehicle モード） |
| `--filter` | 主要道路のみ（国道・都道府県道・高速 or 幅員 5.5m 以上）に絞り込み |
| `--nationwide` | 全国版（`input/` 以下を全件処理） |

### `make_access_links.py` 主なオプション

| オプション | 説明 |
|---|---|
| `--meshes 5338,5339,...` | 対象 1 次メッシュコード |
| `--case {name}` | ケース名（`ksj_to_network_csv.py` と同じ値を指定） |
| `--level 6` | アクセスリンクのメッシュレベル（5=250m, 6=125m） |
| `--pref {都道府県名}` | 都道府県でクリップ |
| `--city {市区町村名}` | 市区町村でクリップ |
| `--mode walk` | 徒歩モード |

---

## `timeslice_search.py` オプション一覧

| オプション | デフォルト | 説明 |
|---|---|---|
| `--links` | `data/saitama/` 道路リンク | 道路リンク parquet パス |
| `--nodes` | `data/saitama/` 道路ノード | 道路ノード parquet パス |
| `--access` | `data/saitama/` アクセスリンク | L6 アクセスリンク parquet パス |
| `--orig-lat` | 35.8578 | 始点緯度 |
| `--orig-lon` | 139.6490 | 始点経度 |
| `--dest-lat` | 36.0420 | 終点緯度 |
| `--dest-lon` | 139.4006 | 終点経度 |
| `--orig-name` | `埼玉県庁` | 始点名（出力ファイル名に使用） |
| `--dest-name` | `東松山市役所` | 終点名（出力ファイル名に使用） |
| `--tmax` | `60` | 移動時間実績（分）。カンマ区切りで複数指定可 |
| `--slice` | `10` | 時刻スライス間隔（分） |
| `--out-dir` | `output/` | 出力ディレクトリ |

---

## 出力ファイル

出力先: `output/`（`--out-dir` で変更可）

### `arrival_map_*.parquet` / `.qml`

始点 A・終点 B それぞれからの到達時間マップ（全 L6 メッシュ）。

| カラム | 内容 |
|---|---|
| `mesh_code` | L6 メッシュコード（11 桁） |
| `到達時間_min` | 最短到達時間（分）。到達不能は NaN |
| `到達時間ランク` | 10 分刻みランク（0〜9）。到達不能は空文字 |
| `geometry` | メッシュポリゴン（EPSG:4326） |

### `timeslice_*.parquet` / `.qml`

通過可能エリア（T_max ごとに 1 ファイル）。

| カラム | 内容 |
|---|---|
| `mesh_code` | L6 メッシュコード（11 桁） |
| `dist_a_min` | 始点 A からの最短到達時間（分）。通過不可は NaN |
| `dist_b_min` | 終点 B への最短到達時間（分）。通過不可は NaN |
| `t_arrive_earliest` | 最早到着可能時刻（= dist_a_min）（分） |
| `t_arrive_latest` | 最遅到着可能時刻（= T_max − dist_b_min）（分） |
| `slack_min` | 余裕時間 = T_max − dist_a − dist_b（分）。通過不可は NaN |
| `slack_rank` | slack の 10 分刻みランク（文字列）。通過不可は空文字 |
| `in_t10` 〜 `in_tXX` | 各時刻スライスで存在可能か（True/False） |
| `geometry` | メッシュポリゴン（EPSG:4326） |

---

## QGIS での可視化

`.qml` ファイルをレイヤーに適用することで色分け表示できる。

### timeslice レイヤーの時刻フィルタ

QGIS の「フィーチャのフィルタ」で以下を指定すると、特定時刻に存在可能なエリアのみ表示できる。

| フィルタ式 | 表示内容 |
|---|---|
| `in_t10 = True` | 出発 10 分後に存在可能なエリア |
| `in_t20 = True` | 出発 20 分後に存在可能なエリア |
| `in_t30 = True` | 出発 30 分後（中間・最も絞り込まれる） |
| `in_t40 = True` | 出発 40 分後に存在可能なエリア |
| `in_t50 = True` | 出発 50 分後に存在可能なエリア |

### カラースキーム

- **arrival_map**: 0〜90 分を 10 分刻みで赤 → 橙 → 黄 → 緑 → シアン の 10 色
- **timeslice**: 通過可能 = 赤、通過不可 = グレーの 2 色

---

## データについて

本ツールは国土交通省「国土数値情報」道路データ（N13）を使用する。

- **出典**: 国土数値情報 道路データ / 国土交通省
- **ダウンロード**: https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N13-v2_1.html
- **ライセンス**: 国土数値情報利用規約に準ずる
- **一方通行**: 国土数値情報には一方通行フィールドが存在しないため、全道路を双方向リンクとして扱っている

### 複製に関する承認

国土数値情報 道路データの原典資料は数値地図（国土基本情報）であり、測量法に基づく国土地理院長承認（複製）**R 6JHf503** を受けている。

> 本製品を複製する場合には、国土地理院の長の承認を得なければならない。
