# 真实 PDB 测试数据来源

下载日期：2026-10-22（America/Los_Angeles）。以下是从 RCSB PDB 下载的完整原始文件，未裁剪、未修改；测试完全离线运行。

| 文件 | 来源与引用 | 下载地址 | SHA-256 |
|---|---|---|---|
| 1CRN.pdb | [RCSB 1CRN](https://www.rcsb.org/structure/1CRN)，crambin。Teeter (1984), PNAS，https://doi.org/10.1073/pnas.81.19.6014 | https://files.rcsb.org/download/1CRN.pdb | `42199a30a0701864a2a5cc76cd7f35cc544cd0e65fbcf63e03c166543249b811` |
| 1IGT.pdb | [RCSB 1IGT](https://www.rcsb.org/structure/1IGT)，免疫球蛋白。Harris et al. (1997), Biochemistry，https://doi.org/10.1021/bi962514+ | https://files.rcsb.org/download/1IGT.pdb | `0a71448e9cb985a3562a4a086f807f86deb235a4f04bdccabf483d5a860b0c95` |

结构 DOI： https://doi.org/10.2210/pdb1CRN/pdb 和 https://doi.org/10.2210/pdb1IGT/pdb 。

选择理由：1CRN 检查普通单链全长与内部裁剪；1IGT 检查多链、真实编号间隔、插入码，以及 HETATM 排除警告。1IGT 的 B/D 链含 52A、82A/B/C、100H/I/J/K。

这些文件是结构坐标验证数据，不是下游设计流程的性能基准。来源、下载版本固定后用哈希防止测试数据无意变化。
