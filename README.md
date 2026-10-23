# Residue Map

A small, local tool for checking protein residue labels, cropping one or more segments of a PDB chain, and exporting an explicit original-to-new residue map.

**Status:** v0.2 prototype. No external Python packages, GPU, account, or API key required. Python 3.10+.

## 快速开始

下载并解压项目。在 VS Code 中打开 `residue-map` 文件夹，再打开 Terminal：

```bash
python3 server.py
```

Windows 如果没有 `python3`，运行 `py server.py`，或使用当前环境的 `python server.py`。
浏览器会自动打开 <http://127.0.0.1:8765>。若未打开，请手动访问该地址。
终端按 Ctrl+C 停止。程序监听本机地址；不要作为公网服务部署。

## 第一次试用

1. 点击 **Try numbering demo**。这个文件是合成编号测试数据，仅有 CA 原子，不是实验结构，不能用于真实蛋白设计。
2. 选择链 A。看到的序列为 `GASK`，原始编号为 `123, 125, 125A, 150`。
3. 在热点框输入 `125A,150`，点击 **Check positions**，应均显示 Found。
4. 选择起点 `125A · SER`、终点 `150 · LYS`，点击 **Prepare downloads**。
5. 映射应为 `125A → 1`、`150 → 2`，链名仍为 A。
6. 下载 `cropped.pdb`、`mapping.csv`、`report.json` 和 `cropped.fasta`。

也可以输入 `125B`，观察工具拒绝不存在的热点。有效但落在裁剪区外的热点会在结果中标为 **outside crop**，不会静默映射到其他位置。

## 新功能：多段裁剪、编号选择与 FASTA

- **多段裁剪**：先选第一段起点/终点，点击 **+ Add range** 添加更多片段；**Remove range** 删除额外片段。所有端点都使用所选链的原始 PDB 标签，区间按实际观察顺序包含两端。
- 不论添加顺序如何，输出都保持原文件顺序；重叠、重复及相邻区间合并，残基不会重复导出。切换链后清除额外区间。
- **Renumber consecutively**：从 **New starting number** 指定的整数开始编号，默认 1，支持 0 和负数；所有生成编号必须位于标准 PDB 的 -999..9999 内。清除插入码，原子序号始终重新生成。
- **Keep original IDs**：保留原始残基编号及插入码（如 125A），起始编号输入停用；仍然生成完整映射。
- **FASTA**：导出所选 ATOM 残基的观察序列，未知名称为 X，每行最多 60 个字母。每个不连续的观察片段分别生成一条 FASTA 记录；不拼接成一条假定连续的序列，不补缺失残基。
- 多段之间在 PDB 中写入 TER；顺序重编号跨片段继续计数。TER 与 FASTA 分段按省略的观察残基判断，不推断真实肽键连续性。

**练习**：演示 A 链，第一段选 `123..123`，添加第二段选 `125A..150`。从 10 编号得到 `123→10、125A→11、150→12`；选择保留原编号则三者编号保持不变。FASTA 有两条记录，序列分别为 `G` 和 `SK`。热点 `125` 标为 outside crop。

## Scope and explicit behavior

- Input: one PDB coordinate model, up to 8 MB in the UI. mmCIF is not supported.
- Only `ATOM` residues are selected. `HETATM` is excluded, **including modified amino acids stored as HETATM**, ligands, ions, and waters. Do not use the exported file when these components must be retained.
- A residue is identified by `(chain ID, residue number, insertion code)`. Blank chain IDs and negative residue numbers are supported.
- Hotspots use original PDB labels, comma-separated, within the selected chain: `125,125A,150`. Ranges are not accepted. They are not sequence-array indices.
- Sequence is derived from observed ATOM records, not SEQRES or a reference sequence. Unknown residue names display as `X`.
- Crop endpoints follow **observed file order**, inclusive. They do not imply consecutive original numbering or continuous peptide connectivity.
- Choose sequential renumbering (default start 1, insertion codes cleared) or preserve original numbers/insertion codes. Chain ID, residue names, coordinates, alternate locations, occupancy and B-factor fields are retained. Atom serials are regenerated.
- Export omits headers, CONECT, SSBOND, ANISOU, HETATM and all other metadata. It does not repair missing atoms/residues, infer bonds or select alternate conformers.
- Multiple models, non-contiguous duplicate residue identifiers and multiple residue names at the same position are rejected rather than guessed.
- Input is sent only to the Python process on your computer. It is kept in memory, not written to disk by the app. Exports are saved when you download them.
- Coordinates and residue mapping are the scope of this prototype. No binding, stability, assay performance or downstream-model compatibility is predicted or guaranteed.

## Downloads

| File | Contents |
|---|---|
| `cropped.pdb` | Selected ATOM records with the chosen numbering; TER separates disjoint observed segments; END terminates the file |
| `mapping.csv` | Original chain, number and insertion code → new identifier |
| `report.json` | Input SHA-256, numbering mode, normalized observed segments, warnings, hotspot mapping and full residue map |
| `cropped.fasta` | Observed ATOM-derived sequence; one record per disjoint selected segment |

`original_observed_index` is a 1-based position among observed residues in the original selected chain. It is **not** a UniProt/SEQRES sequence position.

## Tests

From the project directory:

```bash
python3 -m unittest discover -s tests -v
```

Tests cover insertion codes, gaps, negative IDs, blank chains, crop bounds, alternate locations, HETATM exclusion, malformed coordinates and ambiguity rejection. Real PDB round-trip regressions use locally bundled 1CRN and 1IGT fixtures; source attribution and hashes are in `tests/fixtures/SOURCES.md`. HTTP and actual browser verification results are recorded in `VALIDATION.md`.

Optional front-end async-state regression (Node.js, development only):

```bash
node tests/test_ui_state.cjs
```

Changing the structure, chain, ranges, numbering options or hotspots invalidates previous downloads; delayed responses from an older selection are discarded.

## API compatibility

`POST /api/export` keeps the original `start`/`end` fields as a single-range fallback. New clients may provide `ranges: [[start, end], ...]` using zero-based observed positions, `numbering: "sequential" | "original"`, and integer `start_number` (default 1). Explicit `ranges` takes priority over legacy bounds; invalid or empty ranges return HTTP 400. Defaults retain the previous crop behavior.

The response adds `fasta`. `report.format_version` is now **2**: it adds numbering mode, normalized observed segments, FASTA scope, and hotspot `new_icode`/`new_label`. Existing numeric `new` is retained, but consumers must use `new_label` or `(new, new_icode)` to distinguish preserved insertion codes; outside-crop entries use null. Mapping CSV columns are unchanged.

## Project layout

```text
core.py        PDB parsing, validation, cropping and mapping
server.py      Local HTTP server and API
index.html     Browser interface, no external assets
examples/      Clearly labeled synthetic demo
tests/        Core regression tests
```

## Next useful steps

- Ask 3 protein-design users to try their own structures and report failures.
- Add mmCIF support using a mature structural-biology library.
- Add a 3D residue picker after the mapping behavior is validated.
- Add a named design-tool exporter only after verifying against that tool's documented input format.

This prototype is not affiliated with Rosetta Commons, BindCraft or any protein-design model.
