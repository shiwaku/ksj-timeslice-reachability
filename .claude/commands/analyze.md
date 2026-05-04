# /analyze

時刻スライス型到達圏分析を実行する。

## サンプルデータ（埼玉県）で実行

```bash
python3 src/timeslice_search.py \
  --orig-lat {ORIG_LAT} \
  --orig-lon {ORIG_LON} \
  --orig-name {ORIG_NAME} \
  --dest-lat {DEST_LAT} \
  --dest-lon {DEST_LON} \
  --dest-name {DEST_NAME} \
  --tmax {TMAX}
```

`--links`・`--nodes`・`--access` を省略すると `network/saitama/` のサンプルデータを使用する。

## 独自ネットワークデータで実行

`/generate-network` でネットワーク生成済みの場合:

```bash
python3 src/timeslice_search.py \
  --links network/{CASE}/KSJ_N13-24_{CASE}_道路リンク.parquet \
  --nodes network/{CASE}/KSJ_N13-24_{CASE}_道路ノード.parquet \
  --access network/{CASE}/KSJ_N13-24_{CASE}_アクセスリンク_L6.parquet \
  --orig-lat {ORIG_LAT} \
  --orig-lon {ORIG_LON} \
  --orig-name {ORIG_NAME} \
  --dest-lat {DEST_LAT} \
  --dest-lon {DEST_LON} \
  --dest-name {DEST_NAME} \
  --tmax {TMAX}
```

## 変数の置き換え

- `{ORIG_LAT}` / `{ORIG_LON}`: 始点の緯度・経度（10進数）
- `{ORIG_NAME}`: 始点名（出力ファイル名に使用。例: `横浜駅`）
- `{DEST_LAT}` / `{DEST_LON}`: 終点の緯度・経度
- `{DEST_NAME}`: 終点名
- `{TMAX}`: 移動時間実績（分）。カンマ区切りで複数指定可（例: `60,65,70`）
- `{CASE}`: `/generate-network` で指定したケース名

## 出力先

`output/` に以下が生成される:

- `arrival_map_{ORIG_NAME}.parquet / .qml` — 始点からの到達時間マップ
- `arrival_map_{DEST_NAME}.parquet / .qml` — 終点からの到達時間マップ
- `timeslice_{ORIG_NAME}_{DEST_NAME}_tmax{TMAX}.parquet / .qml` — 通過可能エリア

QGIS で `.qml` をレイヤーに適用し、`in_t30 = True` などのフィルタで特定時刻の存在可能エリアを表示できる。

## 注意

- 始点・終点が対象ネットワーク外の場合はエラーになる
- `--tmax` が最短経路時間より小さい場合も通過可能メッシュが 0 になる
- **一方通行は考慮されていない**（国土数値情報に一方通行フィールドなし）
