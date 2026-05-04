# /generate-network

指定された都道府県・市区町村の道路ネットワークデータを生成する。

## 手順

1. 対象エリアの 1 次メッシュコードを確認する（`README.md` の主要都道府県コード表を参照）
2. `input/` に国土数値情報 GeoJSON が配置済みか確認する
3. 以下のコマンドを実行する

### 都道府県単位

```bash
python3 src/ksj_to_network_csv.py \
  --meshes {MESHES} \
  --case {CASE} \
  --pref {PREF}

python3 src/make_access_links.py \
  --meshes {MESHES} \
  --case {CASE} \
  --level 6 \
  --pref {PREF}
```

### 市区町村単位

```bash
python3 src/ksj_to_network_csv.py \
  --meshes {MESHES} \
  --case {CASE} \
  --city {CITY} \
  --pref {PREF}

python3 src/make_access_links.py \
  --meshes {MESHES} \
  --case {CASE} \
  --level 6 \
  --city {CITY} \
  --pref {PREF}
```

## 変数の置き換え

- `{MESHES}`: 1 次メッシュコード（カンマ区切り。例: `5238,5239,5338,5339`）
- `{CASE}`: ケース名（例: `kanagawa_pref`、`yokohama_city`）
- `{PREF}`: 都道府県名（例: `神奈川県`）
- `{CITY}`: 市区町村名（例: `横浜市`）。`--city` 指定時は `--pref` も必須

## 出力先

`network/{CASE}/` に以下が生成される:

- `KSJ_N13-24_{CASE}_道路リンク.parquet`
- `KSJ_N13-24_{CASE}_道路ノード.parquet`
- `KSJ_N13-24_{CASE}_アクセスリンク_L6.parquet`

生成後は `/analyze` コマンドで分析を実行できる。
