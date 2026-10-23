# Residue Map

**检查热点 · 裁剪结构 · 保留编号映射**  
**Check positions · Crop structures · Keep residue identities traceable**

A lightweight local tool for inspecting PDB residue identities, cropping one or more segments of a protein chain, and exporting explicit original-to-new numbering maps.

轻量级本地 PDB 工具：核对残基位置，裁剪蛋白链的一个或多个片段，并导出可追溯的原编号到新编号映射。

[中文说明](#中文说明) · [English guide](#english-guide)

**v0.2 prototype / 原型版本** · Python 3.10+ · MIT License  
无需第三方 Python 包、GPU、账号或 API Key。  
No third-party Python packages, GPU, account, or API key required.

---

## 中文说明

### 这个工具解决什么问题？

结构中的残基编号可能存在间隔、负数和插入码。裁剪或重新编号后，仅看新文件容易失去原始位置对应关系。

Residue Map 将“检查位置 → 选择片段 → 编号 → 导出映射”放在一个本地页面中。例如：

| 原始链 | 原始编号 | 裁剪后编号 |
|---|---|---|
| A | 125A | 1 |
| A | 150 | 2 |

这里的“热点检查”只确认指定位置是否存在，不预测它是否影响结合或蛋白功能。

### 快速开始

安装 Python 3.10 或更新版本，然后下载本仓库：

```bash
git clone https://github.com/alexzhang0713-rgb/residue-map.git
cd residue-map
python3 server.py
```

也可以在仓库页面点击 **Code → Download ZIP**，解压后在包含 `server.py` 的文件夹打开终端，运行 `python3 server.py`。

Windows 可使用：

```powershell
py server.py
```

浏览器会自动打开 [http://127.0.0.1:8765](http://127.0.0.1:8765)。若未打开，请手动访问；终端按 **Ctrl+C** 停止服务。程序仅监听本机地址，不应直接作为公网服务部署。

可选启动方式：

```bash
python3 server.py --port 8766 --no-browser
```

### 第一次试用

1. 点击 **Try numbering demo**。这是仅含 CA 原子的合成编号数据，不是实验结构，不能用于真实蛋白设计。
2. 选择 A 链，序列应为 `GASK`，原始编号为 `123, 125, 125A, 150`。
3. 输入热点 `125A,150`，点击 **Check positions**，应均显示 **Found**。
4. 起点选择 `125A · SER`，终点选择 `150 · LYS`；保持连续编号、起始编号为 1。
5. 点击 **Prepare downloads**，核对 `125A → 1`、`150 → 2`。
6. 下载 PDB、CSV、JSON 和 FASTA 四种文件。

输入 `125B` 可观察不存在位置的报错。存在但未被选入裁剪范围的热点会标为 **outside crop**。

### 功能与操作

| 功能 | 操作与行为 |
|---|---|
| 结构导入 | 上传单模型 `.pdb` / `.ent` 文件，页面限制为 8 MB |
| 链与序列查看 | 选择一条链，查看文件中实际有 ATOM 坐标的残基序列 |
| 热点核对 | 用逗号分隔原始标签，如 `125,125A,150`；热点输入不支持范围 |
| 多段裁剪 | 用 **+ Add range** 添加区间，用 **Remove range** 删除额外区间 |
| 连续编号 | **Renumber consecutively** 从指定整数开始，默认 1；支持 0 和负数，全部结果须落在 `-999..9999` |
| 原编号保留 | **Keep original IDs** 保留原残基编号与插入码，如 `125A` |
| FASTA 导出 | 每个不连续的观察片段分别生成一条序列记录，未知残基为 `X` |

所有裁剪端点均包含在内，按原文件观察顺序判断。区间按原顺序合并，重叠、重复及相邻区间去重合并。切换链会清除额外区间。

**多段练习**：演示 A 链选择 `123..123` 和 `125A..150`。从 10 编号时得到 `123→10、125A→11、150→12`；保留原编号时标签不变。FASTA 分别包含 `G` 和 `SK`，热点 `125` 标为 outside crop。

不连续片段之间写入 `TER`，连续编号跨片段继续计数。修改结构、链、区间、编号或热点后，旧下载会失效；过期请求不会覆盖新选择。

### 下载文件

| 文件 | 内容 |
|---|---|
| `cropped.pdb` | 所选 ATOM 记录与指定编号，片段间以 TER 分隔，末尾为 END |
| `mapping.csv` | 原链、原编号、插入码与新身份的完整对应关系 |
| `report.json` | 输入 SHA-256、编号模式、观察片段、警告、热点结果及映射 |
| `cropped.fasta` | 所选 ATOM 残基的观察序列，每个不连续片段一条记录，每行最多 60 个字母 |

`original_observed_index` 是原始所选链中**从 1 开始的观察位置**，不是 UniProt 或 SEQRES 位置。建议同时保存原文件和导出文件，以便追踪。

### 适用范围与限制

- 残基身份由 **链名、整数编号、插入码** 共同确定；支持空白链名和负编号。
- 只选择 `ATOM`。排除 `HETATM`，包括配体、离子、水和以 HETATM 保存的修饰氨基酸；需要保留这些组分时，当前导出不适用。
- 只支持单模型 PDB，不支持 mmCIF。多模型、非连续重复身份和同位置多个残基名会被拒绝。
- 序列来自实际观察到的 ATOM，不来自完整参考序列。不会补全缺失残基或原子；未知残基名显示为 `X`。
- 坐标、残基名称、链名、备用原子位置、occupancy 和 B-factor 保留；原子序号重新生成。连续编号清除插入码，原编号模式保留插入码。
- 不保留头部、CONECT、SSBOND、ANISOU 等元数据，不修复结构、不推断键、不自动选择构象。
- 观察顺序、编号间隔及导出分段均不能证明真实肽键连续性。FASTA 也不代表完整蛋白序列。
- 浏览器仅将输入发送至你电脑上的 Python 服务；应用在内存处理，不主动保存输入文件。下载时才保存导出文件。
- 不预测结合、稳定性或实验表现，也不保证特定下游设计工具兼容。

### 测试与开发

```bash
python3 -m unittest discover -s tests -v
```

当前有 27 项 Python 测试，覆盖编号、插入码、空白链、备用原子位置、多段裁剪、编号模式、FASTA 和 HTTP 接口。真实结构回归使用随仓库提供的 1CRN / 1IGT；来源、文献与哈希见 [SOURCES.md](tests/fixtures/SOURCES.md)，测试及浏览器验证记录见 [VALIDATION.md](VALIDATION.md)。

可选前端状态回归需要 Node.js，仅用于开发：

```bash
node tests/test_ui_state.cjs
```

与 pdb-tools 的可复现对比见 [compare_pdbtools.py](tools/compare_pdbtools.py)。pdb-tools 是该开发比较的可选依赖，应用运行不需要它。

---

## English guide

### What does it do?

PDB residue numbering can contain gaps, negative numbers, and insertion codes. Cropping or renumbering a structure can obscure how its new positions relate to the original file.

Residue Map brings position checking, segment selection, numbering, and mapping exports into one local interface. For example, original chain A residues `125A` and `150` can become `1` and `2`, with the original identities recorded in a mapping table.

“Hotspot checking” verifies that a supplied position exists. It does not identify functional hotspots or predict binding effects.

### Quick start

Install Python 3.10 or later, then clone the repository:

```bash
git clone https://github.com/alexzhang0713-rgb/residue-map.git
cd residue-map
python3 server.py
```

Alternatively, choose **Code → Download ZIP**, extract it, and run `python3 server.py` from the directory containing `server.py`.

On Windows:

```powershell
py server.py
```

The app opens [http://127.0.0.1:8765](http://127.0.0.1:8765) automatically. Open it manually if needed. Press **Ctrl+C** in the terminal to stop the server. It listens on localhost; do not deploy it directly as a public service.

To use another port without opening a browser:

```bash
python3 server.py --port 8766 --no-browser
```

### Try the demo

1. Click **Try numbering demo**. This is synthetic CA-only numbering data, not an experimental structure or a valid input for real protein design.
2. Select chain A: sequence `GASK`, original labels `123, 125, 125A, 150`.
3. Enter `125A,150` and click **Check positions**. Both should be **Found**.
4. Select `125A · SER` through `150 · LYS`, with consecutive numbering starting at 1.
5. Click **Prepare downloads** and verify `125A → 1` and `150 → 2`.
6. Download the PDB, CSV, JSON, and FASTA files.

Try `125B` to see a missing-position error. Existing hotspots outside the selection are explicitly marked **outside crop**.

### Features and behavior

| Feature | How it works |
|---|---|
| Import | Upload a single-model `.pdb` / `.ent` file, up to 8 MB in the UI |
| Chain inspection | View the observed ATOM-derived sequence of one selected chain |
| Position checking | Enter comma-separated original labels such as `125,125A,150`; hotspot ranges are not accepted |
| Multiple segments | Use **+ Add range** and **Remove range** to manage additional inclusive ranges |
| Consecutive numbering | Start from an integer of your choice, default 1; zero and negative starts are allowed if all output IDs fit `-999..9999` |
| Original numbering | **Keep original IDs** preserves residue numbers and insertion codes, including labels such as `125A` |
| FASTA export | One observed-sequence record per disjoint selected segment; unknown names become `X` |

Endpoints follow observed file order. Output remains in that order, regardless of range-entry order. Overlapping, duplicate, and adjacent ranges are merged without duplicating residues. Changing the chain clears additional ranges.

**Multiple-segment example:** in demo chain A, select `123..123` and `125A..150`. Starting at 10 yields `123→10`, `125A→11`, `150→12`. Original-numbering mode preserves the labels. FASTA contains separate `G` and `SK` records; hotspot `125` is outside the crop.

Disjoint observed segments are separated by `TER` records. Consecutive numbering continues across segments. Editing the structure, chain, ranges, numbering, or hotspots invalidates previous downloads; stale responses are discarded.

### Downloads

| File | Contents |
|---|---|
| `cropped.pdb` | Selected ATOM records with the chosen numbering, TER between segments, and END |
| `mapping.csv` | Original chain, residue number, and insertion code mapped to the new identity |
| `report.json` | Input SHA-256, numbering mode, observed segments, warnings, hotspots, and full mapping |
| `cropped.fasta` | ATOM-derived observed sequence, one record per disjoint segment, wrapped at 60 letters |

`original_observed_index` is a **1-based observed position in the original selected chain**, not a UniProt or SEQRES position. Keep the original input together with its exports for traceability.

### Scope and limitations

- Residue identity is **chain ID + integer residue number + insertion code**. Blank chain IDs and negative numbers are supported.
- Only `ATOM` residues are selected. `HETATM` is excluded, including ligands, ions, waters, and modified amino acids stored as HETATM. These exports are unsuitable when those components must be retained.
- Input is single-model PDB; mmCIF is not supported. Multiple models, non-contiguous duplicate identities, and multiple residue names at one position are rejected.
- Sequences reflect observed ATOM records, not a complete reference sequence. Missing residues and atoms are not reconstructed; unknown residue names appear as `X`.
- Coordinates, residue names, chain IDs, alternate locations, occupancy, and B-factors are retained. Atom serials are regenerated. Consecutive numbering clears insertion codes; original-numbering mode preserves them.
- Headers, CONECT, SSBOND, ANISOU, and other metadata are omitted. The app does not repair structures, infer bonds, or select alternate conformers.
- Observed order, numbering gaps, and exported segments do not establish peptide connectivity. FASTA is not necessarily a full protein sequence.
- Input is sent only to the Python process on your computer and processed in memory. The app does not write input files to disk; exports are saved when downloaded.
- Binding, stability, experimental performance, and compatibility with specific downstream design tools are not predicted or guaranteed.

### Tests and development

```bash
python3 -m unittest discover -s tests -v
```

The current suite has 27 Python tests covering residue identities, insertion codes, blank chains, alternate locations, multiple ranges, numbering modes, FASTA, and HTTP behavior. Bundled real-structure fixtures are 1CRN and 1IGT. Attribution, literature references, and hashes are in [SOURCES.md](tests/fixtures/SOURCES.md); validation and browser-check records are in [VALIDATION.md](VALIDATION.md).

Optional front-end state regression, requiring Node.js for development only:

```bash
node tests/test_ui_state.cjs
```

A reproducible comparison with pdb-tools is available in [compare_pdbtools.py](tools/compare_pdbtools.py). That comparison has an optional pdb-tools dependency; the app itself does not.

---

## API / 接口

**中文：** `POST /api/export` 支持下列参数。`ranges` 使用链内从 0 开始的观察位置，并优先于旧的 `start` / `end`；旧参数继续支持单区间。空区间或非法参数返回 HTTP 400。

**English:** `POST /api/export` accepts the fields below. `ranges` uses zero-based observed chain positions and takes priority over legacy `start` / `end`, which remain available for a single range. Empty or invalid selections return HTTP 400.

```json
{
  "pdb": "<PDB text>",
  "chain": "A",
  "ranges": [[0, 0], [2, 3]],
  "numbering": "sequential",
  "start_number": 10,
  "hotspots": "123,125A"
}
```

**中文：** `numbering` 为 `sequential` 或 `original`，默认连续编号、从 1 开始。响应包含 `pdb`、`csv`、`fasta`、`mapping` 和 `report`。`report.format_version` 为 **2**，包含编号模式、片段及 FASTA 范围说明。热点保留数值 `new`，完整编号须使用 `new_label` 或 `(new, new_icode)`；裁剪区外为 null。CSV 列保持不变。

**English:** `numbering` is `sequential` or `original`, defaulting to consecutive numbering from 1. The response contains `pdb`, `csv`, `fasta`, `mapping`, and `report`. Report format version **2** includes numbering, segments, and FASTA scope. Hotspots retain numeric `new`; use `new_label` or `(new, new_icode)` for complete identities. Outside-crop values are null. CSV columns are unchanged.

## Project layout / 项目结构

| Path / 路径 | Purpose / 用途 |
|---|---|
| `core.py` | Parsing, selection, numbering, mapping, FASTA / 解析、选择、编号、映射与 FASTA |
| `server.py` | Local HTTP server and API / 本地服务与接口 |
| `index.html` | Browser UI, no external assets / 页面界面，无外部资源 |
| `examples/` | Clearly labeled synthetic demo / 明确标注的合成演示 |
| `tests/` | Core, HTTP, UI-state, and real-structure regressions / 核心、接口、页面状态与真实结构回归 |
| `tools/` | Optional development comparison script / 可选开发对比脚本 |

## Feedback and roadmap / 反馈与后续方向

遇到问题请在 [Issues](https://github.com/alexzhang0713-rgb/residue-map/issues) 提供复现步骤、预期与实际结果，以及可公开分享的最小示例。后续方向包括更多真实用户验证、基于成熟库的 mmCIF 支持、3D 残基选择，以及经过验证的下游工具导出。

Please report reproducible problems through [Issues](https://github.com/alexzhang0713-rgb/residue-map/issues), with expected and actual behavior and a minimal example you can share publicly. Possible next steps include more real-user validation, mmCIF support through a mature library, 3D residue selection, and validated downstream-tool exports.

## License / 许可证

项目代码采用 [MIT License](LICENSE)。真实结构测试数据的来源及引用见 [SOURCES.md](tests/fixtures/SOURCES.md)。

Project code is released under the [MIT License](LICENSE). Attribution and references for real-structure fixtures are listed in [SOURCES.md](tests/fixtures/SOURCES.md).

本项目与 Rosetta Commons、BindCraft 或任何蛋白设计模型无隶属关系。  
This project is not affiliated with Rosetta Commons, BindCraft, or any protein-design model.
