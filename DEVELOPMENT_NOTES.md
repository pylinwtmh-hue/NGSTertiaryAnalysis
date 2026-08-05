# 三級分析 Pipeline 開發筆記（v3.6）

**負責人：** 林伯昱（p88124019@gs.ncku.edu.tw）
**最後更新：** 2026-08-03

> 這份是**內部開發紀錄**（環境、部署、踩雷記錄、驗證過程）。
> 對外的使用說明在 `README.md`，給同事的操作手冊在 `tertiary_pipeline_guide.md`。
>
> 檔名原本是小寫的 `readme.md`，但 Windows 檔案系統不分大小寫，
> 跟 `README.md` 撞名會讓 `git pull` / `git clone` 出錯，
> 2026-08 改名為 `DEVELOPMENT_NOTES.md`。

---

## 目錄

1. [環境規則（必讀）](#環境規則)
2. [目錄結構](#目錄結構)
3. [容器建立](#容器建立)
4. [資料庫下載與建置](#資料庫下載與建置)
5. [Pipeline 執行](#pipeline-執行)
6. [結果驗證](#結果驗證)
7. [PGx Module 建置記錄](#pgx-module-建置記錄)
8. [踩雷記錄（完整版）](#踩雷記錄)
9. [v3.5 更新記錄](#v35-更新記錄)
10. [v3.6 更新記錄](#v36-更新記錄)

---

## 環境規則

```bash
conda activate base      # apptainer build/pull sif 用這個
conda activate nextflow  # 跑 nextflow pipeline 用這個
conda activate genome    # 跑 bcftools / python 工具用這個
```

### Apptainer build 必加參數

mksquashfs 4.7.5 有 bug，`--disable-cache` 和 `APPTAINER_SQUASH_OPTIONS="-processors 1"` **兩個都要加**，缺一不可：

```bash
# build from def file
APPTAINER_SQUASH_OPTIONS="-processors 1" \
apptainer build --disable-cache \
    /data/pylin1991/nf-containers/{sif名稱}.sif \
    /tmp/{def名稱}.def

# pull from Docker Hub
APPTAINER_SQUASH_OPTIONS="-processors 1" \
apptainer pull --disable-cache \
    /data/pylin1991/nf-containers/{sif名稱}.sif \
    docker://image:tag
```

---

## 目錄結構

```
/data/pylin1991/nf-containers/
├── NGStertiary/1_0_0/
│   ├── main_tertiary.nf               ✅ v3.5（依 pipeline_type 建 channel → 組合 sub-workflow）
│   ├── nextflow_tertiary.config       ✅ v3.5
│   ├── modules/
│   │   ├── prepare_vcf.nf             ✅ NCKUH ensemble VCF 前處理
│   │   ├── prepare_vcf_dragen.nf      ✅ DRAGEN VCF 前處理（chrM 分流）+ PLOIDY_REPORT_DRAGEN
│   │   ├── snv_annotation.nf          ✅ VEP 115 + Pangolin
│   │   ├── parse_csq.nf               ✅ transcript 選取 + TSV（61 欄，含 STRAND_BIAS）
│   │   ├── acmg_classifier.nf         ✅ ACMG（ClinGen SVI 2022）
│   │   ├── annotate_snv.nf            ✅ v3.5 SNV 尾段 sub-workflow（SNV_ANNOTATE→PARSE_VEP_CSQ→ACMG）
│   │   ├── mito_annotation.nf         ✅ mtDNA annotation
│   │   ├── str_annotation.nf          ✅ STRchive 分類 + ANNOTATE_STR_{NCKUH,DRAGEN} sub-workflow
│   │   ├── cnv_sv_annotation.nf       ✅ AnnotSV × 5 processes + ANNOTATE_CNV_SV_{NCKUH,DRAGEN} sub-workflow
│   │   └── pgx_annotation.nf          ✅ v3.5（PharmCAT + StellarPGx + OptiType + GATK gVCF + MT-RNR1）
│   ├── scripts/
│   │   ├── add_callers_tag.py
│   │   ├── add_dragen_tag.py
│   │   ├── build_clinvar_lookup.py
│   │   ├── build_gene_moi.py
│   │   ├── build_gnomad_mito_lookup.py
│   │   ├── build_str_lookup.py
│   │   ├── parse_vep_csq.py
│   │   ├── parse_mito_vcf.py
│   │   ├── parse_str_vcf.py
│   │   ├── prepare_sv_dragen.py
│   │   ├── parse_dragen_ploidy.py       ✅ v3.5（DRAGEN ploidy.vcf → ploidy_qc.txt，NDC 與二級統一）
│   │   ├── build_clingen_erepo_lookup.py ✅ v3.6（ERepo → VCEP 判讀查表）
│   │   ├── build_dbnsfp_pknn.py         （dbNSFP + P-KNN 合併；4.9c/5.3a 通用）
│   │   ├── acmg_classifier.py
│   │   ├── parse_pgx_report.py        ✅ v3.4（PharmCAT JSON → pgx.tsv）
│   │   ├── build_outside_calls.py     ✅ v3.4（StellarPGx + OptiType → outside calls）
│   │   ├── parse_stellarpgx.py        ✅ v3.4（StellarPGx .alleles → TSV）
│   │   └── parse_optitype.py          ✅ v3.4（OptiType result.tsv → TSV，Python 3.5 相容）
│   └── stellarpgx_repo/               ✅ git clone SBIMB/StellarPGx
├── tertiary_python_1.0.0.sif
├── vep_115.sif
├── pangolin_cu130_1.0.0.sif            ✅ local 開發機專用（Blackwell sm_120，需驅動 r580+）
├── pangolin_cu121_1.0.0.sif            ✅ production：dgm + dgx（sm_50…sm_90，含 V100 sm_70）
├── annotsv_3.5.10.sif
├── pharmcat_3.2.0.sif                 ✅ v3.3 新增
├── stellarpgx_graphtyper2.5.1.sif     ✅ v3.3 新增
├── samtools_1.23.1.sif                ✅ v3.4 新增（PGx HLA extract 用）
└── optitype_1.3.5.sif                 ✅ v3.4 新增（自建，含 razers3 + samtools）
```

---

## 容器建立

### tertiary_python_1.0.0.sif

```bash
cat > /tmp/tertiary_python.def << 'EOF'
Bootstrap: docker
From: python:3.11-slim

%post
    apt-get update -qq && apt-get install -y --no-install-recommends \
        gcc g++ make \
        zlib1g-dev libbz2-dev liblzma-dev \
        libcurl4-openssl-dev libssl-dev \
        procps bcftools tabix \
        && rm -rf /var/lib/apt/lists/*
    pip install --no-cache-dir cyvcf2 pandas numpy

%test
    python3 -c "import cyvcf2; print('cyvcf2:', cyvcf2.__version__)"
    bcftools --version | head -1
    bgzip --version | head -1

%labels
    Version 1.0.0
    Description "三級分析 Python 工具容器（cyvcf2 + pandas + bcftools + bgzip + tabix）"
EOF

conda activate base
APPTAINER_SQUASH_OPTIONS="-processors 1" \
apptainer build --disable-cache \
    /data/pylin1991/nf-containers/tertiary_python_1.0.0.sif \
    /tmp/tertiary_python.def
```

> ⚠️ `tabix` apt 套件才有 bgzip 和 tabix，只裝 `bcftools` 不夠。

### vep_115.sif

```bash
cat > /tmp/vep_115.def << 'EOF'
Bootstrap: docker
From: ensemblorg/ensembl-vep:release_115.0

%post
    apt-get update -qq && apt-get install -y --no-install-recommends \
        git procps samtools bcftools \
        && rm -rf /var/lib/apt/lists/*

    mkdir -p /opt/vep/Plugins
    mkdir -p /opt/vep/Plugins/loftee_data

    cd /opt/vep/src/ensembl-vep
    perl INSTALL.pl \
        --AUTO p \
        --PLUGINS dbNSFP,LoFtool \
        --PLUGINSDIR /opt/vep/Plugins \
        --NO_HTSLIB \
        --NO_UPDATE \
        2>&1 | tail -10

    cd /tmp
    git clone --depth 1 --branch grch38 \
        https://github.com/konradjk/loftee.git loftee_grch38
    cp -r /tmp/loftee_grch38/* /opt/vep/Plugins/
    rm -rf /tmp/loftee_grch38

    cpanm --quiet --notest List::MoreUtils

%test
    vep --help 2>&1 | grep "ensembl-vep"
    ps --version
    samtools --version | head -1
    bcftools --version | head -1
    perl -e "use DBD::SQLite; print 'DBD::SQLite OK\n'"
    perl -e "use List::MoreUtils; print 'List::MoreUtils OK\n'"
    ls /opt/vep/Plugins/*.pl | wc -l

%labels
    Version 1.0.5
    VEP_release 115
    Description "VEP 115 + LOFTEE GRCh38 complete + samtools + bcftools + dbNSFP + LoFtool"
EOF

conda activate base
apptainer build /data/pylin1991/nf-containers/vep_115.sif /tmp/vep_115.def
apptainer test /data/pylin1991/nf-containers/vep_115.sif
```

### pangolin_cu130_1.0.0.sif ＋ pangolin_cu121_1.0.0.sif（Version 1.0.6）

**同一份 def 檔，改 `TARGET` 一行建出兩顆。** GPU 世代不同、PyTorch wheel 的 arch list
不相容，沒有 prebuilt wheel 能通吃（實測依據見下方 §「為什麼是兩顆容器」）。

#### 📌 兩顆容器對照表（要查哪台用哪顆，看這裡）

| | `pangolin_cu130_1.0.0.sif` | `pangolin_cu121_1.0.0.sif` |
|---|---|---|
| 角色 | **開發機專用**（開發／測試）| **production**（產出臨床報告）|
| def 的 `TARGET` | `cu130` | `cu121` |
| PyTorch | cu130（版本待第一次 build 後 pin）| `torch==2.5.1` **cu121** |
| pip index | `download.pytorch.org/whl/cu130` | `download.pytorch.org/whl/cu121` |
| 必須有的 arch | `sm_120` | `sm_70` |
| 實測 arch flags | （cu128 曾實測 `sm_75 sm_80 sm_86 sm_90 sm_100 sm_120 compute_120`；cu130 建完請補上）| `sm_50 sm_60 sm_70 sm_75 sm_80 sm_86 sm_90` |
| 用在哪些 profile | **`local`** 只有開發機（RTX PRO 6000 Blackwell，sm_120）<br>＝ 唯一被 GPU 世代逼著分家的一台 | **`dgm`**（RTX 2000 Ada，sm_89）<br>**`dgx`**（DGX-2 Tesla V100 ×6，sm_70）<br>＝ **兩台 production 共用這顆** |
| 需要的驅動 | **r580+**（CUDA 13.0 的硬要求）| CUDA 12.1（DGX-2 只有 12.2，cu121 ≤ 12.2 更保險）|

- 切換機制：`nextflow_tertiary.config` 三個 profile 各自宣告 `params.pangolin_sif`，
  `PANGOLIN_SCORE` 用 `container "${params.sif_dir}/${params.pangolin_sif}"`。
  **只能寫在 profile 內**（全域 params 區塊在 `profiles` 之後，會蓋掉 profile 的值）。
- 檔名直接標 CUDA 變體（`cu121` / `cu130`），因為**差異的本質就是 CUDA 變體**，
  不是機器型號 —— 舊名 `pangolin_v100_*` 在 DGM 也用同一顆之後就名不符實了。
- **開發機用 cu130 是刻意的**：Blackwell 是 CUDA 13 世代的原生目標，理論上 13.x 的
  cuBLAS/cuDNN sm_120 kernel 更貼合這張卡。**只有開發機能這樣選** —— cu130 沒有
  `sm_70`，拿去 DGX-2 必爆。
- 分辨手上是哪一顆：`apptainer exec $SIF head -1 /opt/build_versions.txt`
  → `# TARGET=cu130 required_arch=sm_120` 或 `# TARGET=cu121 required_arch=sm_70`

#### 📋 為什麼 DGM 跟 DGX-2 用同一顆 cu121（評鑑可追溯性）

DGM 那張卡兩顆容器都跑得動（見機器表），所以這純粹是政策選擇。選 cu121 是為了
讓評鑑時能講一句乾淨的話：

> **所有 production 分析（DGM + DGX-2）使用同一顆容器映像。**

若 DGM 跟 DGX-2 用不同 torch 版本，就得額外論證「兩個版本產出的 Pangolin splice 分數
等價」—— 那是要拿資料去證的事（同樣本兩台各跑一次、比對分數），為了省一次 rsync
而增加這種舉證責任並不划算。統一成一顆就完全不必談。

開發機（Blackwell / `sm_120`）**必然是例外** —— cu121 沒有 sm_120，物理上不可能統一。
但這不破壞上面那句聲明，因為開發機不產出臨床報告。完整說法是：

> production 兩台（DGM、DGX-2）使用同一顆 cu121 容器；開發機因 GPU 世代限制
> （Blackwell sm_120 不在 cu121 的 arch list 內）使用 cu130 容器，僅用於開發與測試，
> 不產出臨床報告。

⚠️ **連帶的驗證要求**：既然臨床報告在 production 產出，**驗證樣本也應該在 production
容器上跑（或至少複跑一次）**。只在開發機（cu130）驗證、卻在 DGX-2（cu121）發報告，
嚴格講驗證沒有涵蓋實際的 production stack —— 這正是評鑑會問的問題。

#### 版本沿革（三次嘗試，都留著當紀錄）

| Version | torch | 結果 |
|---------|-------|------|
| 1.0.4 | `pip install torch`（**未 pin**，實際解析成 cu130）| 開發機正常 → 以為沒事；部署 DGX-2 才發現 V100 完全不能用（CUDA 13 已移除 Volta）|
| 1.0.5 | `torch==2.5.1+cu121`（單一容器）| DGX-2 修好了，但**反過來把開發機弄壞**（cu121 沒有 sm_120）|
| 1.0.6（過程）| 兩顆：cu128 / cu121 | 試過用 cu128 一顆通吃，被 `%test` 守門員以「缺 sm_70」擋下 → 確認必須兩顆 |
| **1.0.6（現行）** | `pangolin_cu130_*`（開發機）<br>`pangolin_cu121_*`（DGM + DGX-2）| 檔名直接標 CUDA 變體。開發機回到 cu130（Blackwell 原生世代）；production 兩台統一 cu121，便於評鑑聲明「同一顆容器」|

**已知踩雷：**

| 問題 | 原因 | 解決方式 |
|------|------|---------|
| `wget: not found` | 官方 image 沒有 wget | 改用 `curl -sL url -o file` |
| `git: not found` | 官方 image 沒有 git | `apt-get install git` |
| `_Info.__new__() missing type_code` | pyvcf3 的 `_Info` 比原版 pyvcf 多一個必填參數，Pangolin 尚未更新 | `%post` patch 1：補上 `None` |
| `map(int, ...)` crash on `Y`/`R`/`W` | hg38 部分座標含 IUPAC ambiguity code，`one_hot_encode` 只處理 A/C/G/T/N | `%post` patch 2：`re.sub(r'[^01234]', '0', seq)`（同 N，全零 encoding）|
| Pangolin segfault（CSQ 過長） | WGS 的 CSQ 可達 270 KB，Pangolin parse 時爆掉 | module 內先 `bcftools annotate -x INFO/CSQ` |
| Pangolin segfault（alt/random contig） | gencode DB 沒有 `chr*_alt` / `chr*_random` / `chrUn_*` 的 gene model | module 內 `grep -E '^#\|^chr([0-9]+\|[XYM])\t'` 只留標準染色體 |
| **`RuntimeError: ... driver ... too old (found version 12020)`**（DGX-2） | `pip install torch` 沒有 pin，PyPI 預設 wheel 漂移成 cu130（CUDA 13.0），CUDA 13 已移除 Volta 且需驅動 r580+ | DGX-2 改用 **cu121** 容器（`pangolin_cu121_1.0.0.sif`），見 §「為什麼是兩顆容器」 |
| **`CUDA error: no kernel image is available for execution on the device`** | wheel 缺這張卡的 arch。兩個方向都發生過：cu121 沒有 `sm_120`（開發機爆）、cu128/cu130 沒有 `sm_70`（DGX-2 爆）| 一台一顆容器，由 `params.pangolin_sif` 切換。**沒有** prebuilt wheel 同時含兩者 |

#### ⚠️ 為什麼是兩顆容器（2026-08 定案，含實測依據）

這個 pipeline 有三台目標機器，但**只有兩個 GPU 世代分組** —— DGX-2 的 V100 自己一組，
其餘（開發機、DGM）都是新世代，**差了三個世代**：

| 機器 | GPU | compute capability | 容器 |
|------|-----|-------------------|------|
| DGX-2（production，`-profile dgx`） | Tesla V100 ×6 | **sm_70**（Volta）| `pangolin_cu121_1.0.0.sif` |
| 開發機（`-profile local`） | RTX PRO 6000 Blackwell Max-Q | **sm_120**（Blackwell）| `pangolin_cu130_1.0.0.sif` |
| DGM Server（`-profile dgm`） | RTX 2000 Ada Generation（AD107 / Ada）<br>driver 580.173.02 | **sm_89** ⚠️ 不在 cu121 的 arch list 內，靠同 major 版本相容（見下）| `pangolin_cu121_1.0.0.sif`（刻意與 dgx 一致，見下）|

#### ⚠️ DGM 的 sm_89 不在 cu121 的 arch list 裡（靠同 major 版本相容）

實測（`nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv`）：

```
NVIDIA RTX 2000 Ada Generation, 8.9, 580.173.02
```

`sm_89` **沒有**出現在 cu121 wheel 的 arch flags（`sm_50 sm_60 sm_70 sm_75 sm_80 sm_86
sm_90`）裡。它能跑是因為 **CUDA cubin 在同一個 major compute capability 版本內向前相容**
—— `sm_86` 的 cubin 可以在 `sm_87` / `sm_89` 上執行（都是 major 8）。這是 NVIDIA 有文件
保證的行為，不是 workaround。

對照一下就知道為什麼 V100 和 Blackwell 沒有這個好處：`sm_70`（major 7）與 `sm_120`
（major 12）跟其他世代**跨了 major**，所以完全不相容 —— 這也正是必須分兩顆容器的根因。

| 卡 | compute cap | cu121 wheel 怎麼滿足它 |
|----|-------------|----------------------|
| DGX-2 Tesla V100 | sm_70（major 7）| arch list 直接有 `sm_70` |
| DGM RTX 2000 Ada | sm_89（major 8）| **沒有 `sm_89`，用 `sm_86` cubin 向前相容** |
| 開發機 RTX PRO 6000 | sm_120（major 12）| ❌ 完全沒有，必須另一顆 cu130 |

**所以 DGM 的部署驗證比另外兩台更重要**：另外兩台是 arch list 直接命中，DGM 是靠相容性
推論。推論可能有例外（例如某個 cuDNN kernel 只出特定 arch），所以 **一定要在 DGM 上跑一次
真的 conv1d forward**（見下方部署驗證步驟），不能只靠這段論述。

> 附帶資訊：DGM 的驅動是 **580.173.02**，已達 CUDA 13.0 的門檻（r580+），所以 DGM
> 技術上也跑得動 cu130 容器 —— 這與「DGM 之前跑 cu130 正常」一致。仍然選 cu121 是
> **政策決定**（與 DGX-2 統一），不是能力限制。

##### 那要不要為 sm_89 再包第三顆（cu124 之類的）？→ 不要

直覺會想「找一個 arch list 裡有 `sm_89` 的 wheel」，但**目前實測過的兩顆都沒有 `sm_89`**：

| wheel | 實測 arch flags | 有 sm_89？ |
|-------|----------------|-----------|
| `torch 2.5.1+cu121` | `sm_50 sm_60 sm_70 sm_75 sm_80 sm_86 sm_90` | ❌ |
| `torch 2.7.0+cu128` | `sm_75 sm_80 sm_86 sm_90 sm_100 sm_120 compute_120` | ❌ |

CUDA **11.8 起就支援 sm_89**，所以 cu121／cu124／cu126／cu128 在工具鏈層面全都支援它 ——
但 PyTorch 官方 wheel **刻意不單獨編 `sm_89`**，因為 `sm_86` 的 cubin 本來就能跑 Ada，
多編一份只是讓 wheel 變大。也就是說：**PyTorch 自己就是靠這個相容性在支援 Ada 卡的**，
我們不是在用什麼旁門左道。

所以再包一顆 cu124 只會得到「同樣沒有 sm_89、只是 torch 版本不同」的容器，
卻要付出：production 兩台不再是同一顆容器（正是上一節花力氣統一掉的東西）。
**淨損。**

真的想要 arch list 命中 `sm_89`，唯一辦法是**從源碼編**並指定
`TORCH_CUDA_ARCH_LIST="7.0;7.5;8.0;8.6;8.9;9.0"`。1.5–3 小時 + 30 GB 暫存，
換到的只是「本來就跑得動的東西改成原生命中」，效能差異對 Pangolin 這種小 conv
可以忽略。**不做。**

判斷依據：**先跑 DGM 上那個 conv1d forward。** 過了就代表相容性成立、案子結案；
真的失敗才需要重新考慮（而且那時候正解是源碼編，不是換一顆 prebuilt wheel）。
若哪天想確認某個候選 wheel 到底有沒有 sm_89，兩分鐘就能查，不必先包容器：

```bash
python3 -m venv /tmp/archchk && /tmp/archchk/bin/pip install -q \
    --index-url https://download.pytorch.org/whl/cu124 torch
/tmp/archchk/bin/python -c \
    "import torch; print(torch.__version__, torch.version.cuda); print(torch._C._cuda_getArchFlags())"
rm -rf /tmp/archchk
```


**為什麼不能一顆通吃（已實測確認）。** CUDA 工具鏈層面，12.8/12.9 是唯一同時支援
sm_70（deprecated 但可編譯）與 sm_120（12.8 才有）的版本 —— ≤12.6 沒有 sm_120，
13.0 起移除了 Volta。看起來 cu128 應該可以通吃。

但**「工具鏈支援」≠「wheel 裡有」**：PyTorch 為了縮小 wheel 會自己砍舊 arch。
實際去問 `torch._C._cuda_getArchFlags()`：

| wheel | 實測 arch flags | sm_70 | sm_120 |
|-------|----------------|-------|--------|
| `torch 2.5.1+cu121` | `sm_50 sm_60 sm_70 sm_75 sm_80 sm_86 sm_90` | ✅ | ❌ |
| `torch 2.7.0+cu128` | `sm_75 sm_80 sm_86 sm_90 sm_100 sm_120 compute_120` | ❌ | ✅ |
| `torch (cu130 預設)` | （CUDA 13 已移除 Volta）| ❌ | ✅ |

→ **沒有任何 prebuilt wheel 同時涵蓋兩者。** cu128 在工具鏈上支援 sm_70，但官方
wheel 沒把它編進去。所以一台一顆容器，由 `params.pangolin_sif` 依 profile 切換。

那個 `compute_120` 也順便否證了 PTX JIT 這條路：PTX **只能往新的 arch 前向相容**，
cu128 只帶 `compute_120` 的 PTX，對 sm_70 完全沒用。

> **唯一能通吃的做法是從源碼編**（base 換 `nvidia/cuda:12.8.0-devel-*`，
> `TORCH_CUDA_ARCH_LIST="7.0;7.5;8.0;8.6;9.0;12.0"`）。要 1.5–3 小時、30+ GB 暫存，
> 容器也會肥很多 —— 為了省一顆 sif 不值得，故不採用。

**根因不是「某次改動被改掉」，而是依賴漂移。** 原始 def 寫的是

```bash
pip install --no-cache-dir torch torchvision     # ← 沒有 pin
```

這行一直都沒有 pin。PyPI 的 `torch` 預設 wheel 原本是 CUDA 12.x，所以當初包起來
在開發機上一切正常；等 PyTorch 把預設換成 cu130 之後，**同一份 def 檔重建就產出不同的
容器**。而且 cu130 在開發機（Blackwell）上是**正常的** —— 這就是為什麼本地測試一路綠燈，
只有部署到 V100 才爆。`%labels` 的 `cu130` 是當時觀察到的結果，正好留下漂移的證據。

**為什麼 build 的時候沒被抓到**：原始 `%test` 只是

```bash
python3 -c "import torch; print('CUDA:', torch.cuda.is_available())"
```

`apptainer build` / `apptainer test` 都**不帶 `--nv`**，所以這行在任何機器上都印 `False`
而且 **exit code 是 0**，測試照過。新版 `%test` 改成不需要 GPU 也能驗的硬檢查：
`%post` 把這顆容器的目標 arch 寫進 `/opt/required_arch.txt`，`%test` 讀出來比對
wheel 的 arch flags，不符就讓 build 爆掉。**cu128 通吃的嘗試就是被這個守門員擋下來的**
（訊息直接指名缺 `sm_70`／DGX-2），沒有再一次「部署到某一台才發現」。

> **SIF 檔名對應 `params.pangolin_sif`**，改名要同步改 `nextflow_tertiary.config`
> 三個 profile 的宣告（`local` / `dgm` / `dgx`）。
> **只能寫在 profile 內**，不可放全域 params 區塊 —— 那個區塊在 `profiles` 之後，
> 會把 profile 的值蓋掉（見「global params > profile params」踩雷記錄）。
>
> **不要走這兩條**：升驅動解決不了 arch list 的問題（cu130 升上去只會從 `driver too old`
> 變成 `no kernel image`，且 V100 永久不能用）；PTX JIT（`CUDA_FORCE_PTX_JIT=1`）
> 也不行 —— 理由見上面 `compute_120`，而且 cuDNN/cuBLAS 只出 cubin 不出 PTX，
> conv1d 照樣沒有 kernel。

#### 重建指令

```bash
cat > /tmp/pangolin.def << 'EOF'
Bootstrap: docker
From: python:3.11-slim

%post
    apt-get update -qq && apt-get install -y --no-install-recommends \
        gcc g++ make git \
        zlib1g-dev libbz2-dev liblzma-dev \
        libcurl4-openssl-dev libssl-dev \
        procps bcftools tabix \
        && rm -rf /var/lib/apt/lists/*

    # ══════════════════════════════════════════════════════════════════
    #  ★★ 要換目標只改這一行 ★★（TARGET 名稱＝輸出的 sif 檔名）
    #     cu121 → production：DGM + DGX-2（V100 sm_70）→ pangolin_cu121_1.0.0.sif
    #     cu130 → 開發機 RTX PRO 6000 Blackwell（sm_120）→ pangolin_cu130_1.0.0.sif
    #  ══════════════════════════════════════════════════════════════════
    TARGET=cu130

    # ── ⚠️ torch 的 CUDA 變體一定要用 --index-url 明確指定 ─────────────────
    #   為什麼不能用預設的 `pip install torch`：
    #     PyPI 的預設 wheel 會隨時間漂移（實際觀察到從 CUDA 12.x 變成 cu130）。
    #     這就是原本那顆容器的 bug：同一份 def 檔重建就產出不同的容器，
    #     而 cu130 沒有 sm_70 → DGX-2 的 V100 永遠跑不起來。
    #     ★ 重點不是「不能用 cu130」，而是「不能讓它隱性漂移」。
    #   為什麼要兩顆而不是一顆通吃（實測結論，不是猜的）：
    #     torch 2.5.1+cu121 arch = sm_50 sm_60 sm_70 sm_75 sm_80 sm_86 sm_90
    #     torch 2.7.0+cu128 arch = sm_75 sm_80 sm_86 sm_90 sm_100 sm_120 compute_120
    #     CUDA 12.8 在工具鏈上支援 sm_70，但官方 wheel 沒把它編進去
    #     → 沒有任何 prebuilt wheel 同時有 sm_70 和 sm_120，只能一台一顆。
    #     （通吃只能從源碼編，要 1.5–3 小時 + 30GB 暫存，不值得。）
    #   torchvision 已移除：Pangolin 不 import 它，留著只是多一個版本配對的束縛。
    if [ "$TARGET" = "cu121" ]; then
        # production 兩台共用。cu121 <= DGX-2 驅動的 CUDA 12.2 → 連 minor-version
        # 相容都不必依賴；2.5.1 是還有 cu121 wheel 的最後一個版本
        # （2.6 起改成 cu118/cu124/cu126）。
        pip install --no-cache-dir \
            --index-url https://download.pytorch.org/whl/cu121 \
            torch==2.5.1
        echo "sm_70" > /opt/required_arch.txt
    else
        # 開發機專用。Blackwell 是 CUDA 13 世代的原生目標，cuBLAS/cuDNN 的
        # sm_120 kernel 在 13.x 較新，理論上比 cu128 更貼合這張卡。
        #   ⚠️ 版本沒有 pin 到 patch level：cu130 的 torch 版本組合我沒有實測過，
        #      這裡只 pin CUDA 變體（＝真正會漂移的那一項）。
        #      **第一次 build 完請看 /opt/build_versions.txt 印出的 torch 版本，
        #      把它寫回這裡變成 torch==X.Y.Z**，之後重建才完全可重現。
        #   ⚠️ cu130 需要驅動 r580+，且沒有 sm_70 —— 絕對不能拿去 DGX-2。
        #      真的誤用了，%test 的守門員會擋（TARGET=cu121 時要求 sm_70）。
        pip install --no-cache-dir \
            --index-url https://download.pytorch.org/whl/cu130 \
            torch
        echo "sm_120" > /opt/required_arch.txt
    fi
    echo "$TARGET" > /opt/build_target.txt
    echo "[BUILD] TARGET=$TARGET required_arch=$(cat /opt/required_arch.txt)"

    pip install --no-cache-dir gffutils biopython pandas pyfastx pyvcf3
    pip install --no-cache-dir git+https://github.com/tkzeng/Pangolin.git

    # ── 記錄實際解析到的版本（unpinned 依賴漂移就是這個容器踩過的坑）──────
    #   寫進 image，之後 `apptainer exec $SIF cat /opt/build_versions.txt` 就能查，
    #   評鑑要求的「版本可追溯」也用得上。兩顆容器長得很像，這也是分辨的依據。
    {
        echo "# TARGET=$(cat /opt/build_target.txt) required_arch=$(cat /opt/required_arch.txt)"
        pip freeze | grep -iE "^(torch|pangolin|pyvcf3|gffutils|pyfastx|biopython|pandas)"
    } > /opt/build_versions.txt
    echo "--- build_versions.txt ---"; cat /opt/build_versions.txt

    # ── Patch pangolin.py：補上 pyvcf3 新增的 type_code 參數 ──────
    # pyvcf3 的 _Info.__new__() 比原版 pyvcf 多一個必填參數 type_code，
    # Pangolin 尚未更新，補上 None（等同預設行為）即可正常運作。
    # Patch 2：one_hot_encode 加入 IUPAC code 處理，
    #   reference genome 部分座標含 Y/R/W 等 IUPAC code，
    #   原本只處理 A/C/G/T/N，其餘字元讓 map(int,...) crash，
    #   改用 re.sub 把所有非數字字元換成 '0'（同 N，全零 encoding）
    # 兩個 patch 都在 pattern 找不到時 raise SystemExit(1) → 上游改版會讓 build
    #   立刻失敗，而不是產出一個會在跑 case 時才爆的容器。
    python3 - <<'PYEOF'
path = "/usr/local/lib/python3.11/site-packages/pangolin/pangolin.py"
with open(path, "r") as f:
    content = f.read()

# Patch 1：_Info type_code
old1 = "\"Format: gene|pos:score_change|pos:score_change|warnings,...\",\'.\',\'.\')"
new1 = "\"Format: gene|pos:score_change|pos:score_change|warnings,...\",\'.\',\'.\', None)"
if old1 in content:
    content = content.replace(old1, new1)
    print("Patch 1 (_Info type_code) OK")
else:
    print("ERROR: Patch 1 pattern not found")
    for i, l in enumerate(content.splitlines()[238:252], 239):
        print(f"{i}: {l}")
    raise SystemExit(1)

# Patch 2：one_hot_encode IUPAC ambiguity code
old2 = "    seq = seq.replace('G', '3').replace('T', '4').replace('N', '0')"
new2 = ("    seq = seq.replace('G', '3').replace('T', '4').replace('N', '0')\n"
        "    import re; seq = re.sub(r'[^01234]', '0', seq)  # IUPAC ambiguity -> N")
if old2 in content:
    content = content.replace(old2, new2)
    print("Patch 2 (one_hot_encode IUPAC) OK")
else:
    print("ERROR: Patch 2 pattern not found")
    for i, l in enumerate(content.splitlines()[15:30], 16):
        print(f"{i}: {l}")
    raise SystemExit(1)

with open(path, "w") as f:
    f.write(content)
print("All patches written OK")
PYEOF

%test
    pangolin --help 2>&1 | head -3

    # ── ★ 這個容器最重要的守門員：wheel 有沒有這台機器的 arch ★ ────────────
    #   關鍵設計：這段**不需要 GPU** 也能驗，所以在 build 階段就會爆，而不是
    #   等部署到某一台才發現。原始 def 那行 `print(torch.cuda.is_available())`
    #   在 build 時永遠印 False 且 exit 0 —— cu130 就是這樣混過自己的測試。
    #   要驗哪個 arch 由 %post 寫進 /opt/required_arch.txt（跟著 TARGET 走），
    #   所以改 TARGET 就等於同時改了 pip 來源與這裡的驗收標準，不會忘記對齊。
    #   （注意：`apptainer build --notest` 會跳過這段，重建時不要加。）
    python3 - <<'PYEOF'
import sys, torch

target   = open("/opt/build_target.txt").read().strip()
required = open("/opt/required_arch.txt").read().strip()
MACHINE  = {"sm_70": "production：DGX-2 Tesla V100（DGM 的 sm_89 靠 sm_86 相容）",
            "sm_120": "開發機 RTX PRO 6000 Blackwell"}

cuda = torch.version.cuda or ""
print(f"TARGET={target} required_arch={required}")
print("torch:", torch.__version__, "| built with CUDA:", cuda)

# 讀編譯期就烙進 binary 的 arch flags。
#   不要用 torch.cuda.get_arch_list()：它在 is_available() 為 False 時回傳空 list，
#   在沒有 GPU 的 build 環境裡問不出東西。
flags = ""
try:
    flags = torch._C._cuda_getArchFlags() or ""
except Exception:
    pass
if not flags:
    flags = torch.__config__.show()
print("arch flags:", flags)

# 檢查 1（只對 production／sm_70）：CUDA 13.x 移除了 Volta 且要求驅動 r580+，
#   DGX-2 的驅動只到 12.2 → TARGET=cu121 誤裝成 13.x 要立刻擋掉。
#   TARGET=cu130 不受此檢查（開發機本來就是 CUDA 13 世代）。
if required == "sm_70" and not cuda.startswith("12"):
    sys.exit(f"FATAL: torch built with CUDA {cuda!r}。TARGET=cu121 需要 CUDA 12.x"
             "（13.x 已移除 Volta/sm_70 且要求驅動 r580+）。"
             "檢查 --index-url 是不是誤指到 cu130。")

# 檢查 2：★ 真正的把關 ★ wheel 裡到底有沒有這台機器的 arch。
#   「CUDA 12.8 支援 sm_70」不等於「這顆 wheel 包了 sm_70」——
#   PyTorch 為了縮小 wheel 會自己砍舊 arch，只有這裡問得到實話。
#   PTX（compute_XX）只能往新 arch 前向相容，所以也接受同號的 compute_。
if required not in flags and required.replace("sm_", "compute_") not in flags:
    sys.exit(f"FATAL: 這個 torch build 缺少 {required}"
             f"（{MACHINE.get(required, target)}）。\n"
             f"       實際 arch flags：{flags}\n"
             "       這台機器會噴 'no kernel image is available for execution "
             "on the device'。\n"
             "       確認 TARGET 與 pip 的 --index-url 是否對應；\n"
             "       詳見 DEVELOPMENT_NOTES.md §「為什麼是兩顆容器」。")

print(f"OK: {required} 在 arch flags 裡 -> {MACHINE.get(required, target)} 可用")
PYEOF

    cat /opt/build_versions.txt
    python3 -c "import gffutils; print('gffutils OK')"
    bcftools --version | head -1
    bgzip --version | head -1
    tabix --version | head -1
    ps --version
    python3 - <<'PYEOF'
import vcf, re
# Patch 1：_Info type_code
info = vcf.parser._Info('TEST', '.', 'String', 'test', '.', '.', None)
print("Patch 1 (_Info type_code) OK")
# Patch 2：one_hot_encode IUPAC
import sys
sys.path.insert(0, '/usr/local/lib/python3.11/site-packages')
from pangolin.pangolin import one_hot_encode
import numpy as np
result = one_hot_encode('ACGTYN', '+')  # Y 是 IUPAC，不應 crash
print("Patch 2 (one_hot_encode IUPAC) OK")
PYEOF

%labels
    Version 1.0.6
    Description "Pangolin tkzeng/Pangolin + PyTorch pinned by TARGET (cu121 or cu130; see /opt/build_versions.txt) + pyvcf3 patched + bcftools + tabix/bgzip"
EOF

conda activate base

# ── 兩顆都要建。差別只有 def 檔裡的 TARGET 那一行 ────────────────────────
# ⚠️ 不要加 --notest，%test 的 arch 守門員要在 build 階段生效

SIFDIR=/data/pylin1991/nf-containers

# (A) 開發機專用（TARGET=cu130 / sm_120）—— def 檔預設就是這個
apptainer build $SIFDIR/pangolin_cu130_1.0.0.sif /tmp/pangolin.def
apptainer test  $SIFDIR/pangolin_cu130_1.0.0.sif
# 期望：OK: sm_120 在 arch flags 裡 -> 開發機 RTX PRO 6000 Blackwell 可用
#
# ★ 第一次建完，把實際的 torch 版本 pin 回 def 檔（cu130 那支 pip 沒有 pin 到 patch）
apptainer exec $SIFDIR/pangolin_cu130_1.0.0.sif grep '^torch' /opt/build_versions.txt
#   例如印出 torch==2.9.0+cu130 → 把 def 裡 cu130 分支的 `torch` 改成 `torch==2.9.0`

# (B) production 兩台共用：DGM + DGX-2（TARGET=cu121 / sm_50…sm_90）
sed -i 's/^    TARGET=cu130$/    TARGET=cu121/' /tmp/pangolin.def
grep -n '^    TARGET=' /tmp/pangolin.def          # 確認真的改到了
apptainer build $SIFDIR/pangolin_cu121_1.0.0.sif /tmp/pangolin.def
apptainer test  $SIFDIR/pangolin_cu121_1.0.0.sif
# 期望：OK: sm_70 在 arch flags 裡 -> production：DGM + DGX-2 Tesla V100 可用
#   （守門員只驗最嚴格的 sm_70；DGM 那張卡也在 cu121 的 arch list 內，見機器表）

# 分辨兩顆容器（TARGET 寫在 build_versions.txt 第一行）
apptainer exec $SIFDIR/pangolin_cu130_1.0.0.sif head -1 /opt/build_versions.txt
apptainer exec $SIFDIR/pangolin_cu121_1.0.0.sif head -1 /opt/build_versions.txt
# 期望：# TARGET=cu130 required_arch=sm_120  ／  # TARGET=cu121 required_arch=sm_70
```

#### 部署與驗證

**每顆容器都要在它的目標機器上實測** —— 這次的教訓就是「只在一台上測會漏」：

```bash
# ── 1) 開發機（Blackwell / sm_120）：用 pangolin_cu130_1.0.0.sif（cu130）────
SIF=/data/pylin1991/nf-containers/pangolin_cu130_1.0.0.sif
apptainer exec --nv $SIF python3 -c \
    "import torch; print('cuda', torch.cuda.is_available(), torch.cuda.get_device_name(0))"
# 期望：cuda True NVIDIA RTX PRO 6000 Blackwell Max-Q Workstation Edition

# ── 2) 把 cu121 那顆傳到兩台 production（不要傳 pangolin_cu130_1.0.0.sif 過去）──
#      DGX-2
rsync -avz --progress \
    /data/pylin1991/nf-containers/pangolin_cu121_1.0.0.sif \
    n101569@10.11.33.75:/datalake_Intermediate/pipeline/nextflow_containers/
#      DGM
rsync -avz --progress \
    /data/pylin1991/nf-containers/pangolin_cu121_1.0.0.sif \
    n101569@192.168.84.91:/home/pipeline/nextflow_containers/

# ── 3) DGX-2（V100 / sm_70）：用 pangolin_cu121_1.0.0.sif ───────────────
ssh n101569@10.11.33.75
SIF=/datalake_Intermediate/pipeline/nextflow_containers/pangolin_cu121_1.0.0.sif
apptainer exec --nv $SIF python3 -c \
    "import torch; print('cuda', torch.cuda.is_available(), torch.cuda.get_device_name(0))"
# 期望：cuda True Tesla V100-SXM3-32GB

# ── 3b) DGM（RTX 2000 Ada / sm_89）：同一顆 pangolin_cu121_1.0.0.sif ─────
#      ★ 這台最需要實測：sm_89 不在 cu121 的 arch list 內，是靠 sm_86 cubin
#        同 major 相容才跑得動（推論，不是命中）→ 下面那個 conv1d forward 必跑
ssh n101569@192.168.84.91
SIF=/home/pipeline/nextflow_containers/pangolin_cu121_1.0.0.sif
apptainer exec --nv $SIF python3 -c \
    "import torch; print('cuda', torch.cuda.is_available(), torch.cuda.get_device_name(0))"
# 期望：cuda True NVIDIA RTX 2000 Ada Generation
apptainer exec --nv $SIF python3 -c \
    "import torch; x=torch.randn(1,4,64,device='cuda'); \
     c=torch.nn.Conv1d(4,8,3).cuda(); print('conv OK', c(x).shape)"
# 期望：conv OK torch.Size([1, 8, 62])   ← 這行過了才算 sm_86→sm_89 相容真的成立

# ── 4) 版本追溯（第一行是 TARGET，可確認拿到的是哪一顆）────────────────
#      production 兩台印出來必須完全一樣 —— 這就是統一容器要拿去講的證據
apptainer exec $SIF cat /opt/build_versions.txt
```

> ⚠️ **`torch.cuda.is_available()` 印 True 還不夠。** arch 不符時它照樣回 True，
> 要等第一次跑 kernel（Pangolin 的 `conv1d`）才噴 `no kernel image is available` ——
> 這正是為什麼手動測 DGX-2 是 True、進 pipeline 才死。**每台都要跑一次真的 forward**：
> ```bash
> apptainer exec --nv $SIF python3 -c \
>   "import torch; x=torch.randn(1,4,64,device='cuda'); \
>    c=torch.nn.Conv1d(4,8,3).cuda(); print('conv OK', c(x).shape)"
> ```
> 印出 `conv OK torch.Size([1, 8, 62])` 才算真的通。

**在容器修好之前的臨時解法**：`--use_gpu_pangolin false` 走 CPU。
`PANGOLIN_SCORE` 的 CPU 路徑已完整處理（不加 `--nv`、不呼叫 `gpu_lock.sh`、
`CUDA_VISIBLE_DEVICES=""`），**結果與 GPU 完全相同**，只是慢；且 Pangolin 只跑
`INFO` 含 `splice` 的候選變異（非全部變異），實務上可能可以接受。
若要長期走 CPU，記得補 `dgx` profile 的 `withLabel: process_gpu { cpus = N }`
並把 `maxForks` 從 6 調小 —— 否則 6 個 CPU inference 會互搶核心。

---

### annotsv_3.5.10.sif

```bash
cat > /tmp/annotsv_3.5.10.def << 'EOF'
Bootstrap: docker
From: ubuntu:22.04

%post
    apt-get update && apt-get install -y \
        tcl tcllib tclx \
        wget curl git \
        bedtools vcftools \
        && apt-get clean

    cd /opt
    wget https://github.com/lgmgeo/AnnotSV/archive/refs/tags/v3.5.10.tar.gz
    tar xzf v3.5.10.tar.gz
    cd AnnotSV-3.5.10 && make PREFIX=/usr/local install
    rm /opt/v3.5.10.tar.gz

%test
    AnnotSV --help 2>&1 | head -3

%labels
    Version 1.0.0
    AnnotSV 3.5.10
EOF

APPTAINER_SQUASH_OPTIONS="-processors 1" \
apptainer build --disable-cache \
    /data/pylin1991/nf-containers/annotsv_3.5.10.sif \
    /tmp/annotsv_3.5.10.def
```

### pharmcat_3.2.0.sif

```bash
conda activate base
APPTAINER_SQUASH_OPTIONS="-processors 1" \
apptainer pull --disable-cache \
    /data/pylin1991/nf-containers/pharmcat_3.2.0.sif \
    docker://pgkb/pharmcat:3.2.0
```

### stellarpgx_graphtyper2.5.1.sif

```bash
# ⚠️ GitHub repo 裡的 containers/stellarpgx-dev.sif 是 Git LFS 指標，不是真實容器
# 必須從 Docker Hub pull
APPTAINER_SQUASH_OPTIONS="-processors 1" \
apptainer pull --disable-cache \
    /data/pylin1991/nf-containers/stellarpgx_graphtyper2.5.1.sif \
    docker://sbimb/stellarpgx:graphtyper2.5.1
```

### samtools_1.23.1.sif

```bash
APPTAINER_SQUASH_OPTIONS="-processors 1" \
apptainer pull --disable-cache \
    /data/pylin1991/nf-containers/samtools_1.23.1.sif \
    docker://quay.io/biocontainers/samtools:1.23.1--ha83d96e_0
```

### optitype_1.3.5.sif（自建）

為什麼自建：
- `fred2/optitype:latest`（v1.3.1）→ Python 2.7/3.5，razers3 不支援 fastq 輸出，`--no-home` 問題
- `quay.io/biocontainers/optitype:1.5.0` → Pyomo 6.10 bug，constraint infeasible

```bash
cat > /tmp/optitype_1.3.5.def << 'EOF'
Bootstrap: docker
From: ubuntu:22.04

%labels
    maintainer p88124019@gs.ncku.edu.tw
    version 1.3.5
    description OptiType 1.3.5 HLA typing with razers3 3.5.12 + samtools 1.21

%environment
    export PATH="/opt/conda/bin:$PATH"
    export DEBIAN_FRONTEND=noninteractive

%post
    export DEBIAN_FRONTEND=noninteractive
    export PATH="/opt/conda/bin:$PATH"

    apt-get update && apt-get install -y wget bzip2 ca-certificates \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

    wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
        -O /tmp/miniforge.sh
    bash /tmp/miniforge.sh -b -p /opt/conda
    rm /tmp/miniforge.sh

    /opt/conda/bin/conda config --add channels defaults
    /opt/conda/bin/conda config --add channels bioconda
    /opt/conda/bin/conda config --add channels conda-forge
    /opt/conda/bin/conda config --set channel_priority strict

    /opt/conda/bin/conda install -y \
        python=3.11 \
        optitype=1.3.5 \
        samtools=1.21 \
        coincbc \
        && /opt/conda/bin/conda clean -afy

%test
    OptiTypePipeline.py --help 2>&1 | head -3
    razers3 --version 2>&1 | head -2
    samtools --version | head -1

%runscript
    exec OptiTypePipeline.py "$@"
EOF

APPTAINER_SQUASH_OPTIONS="-processors 1" \
apptainer build --disable-cache \
    /data/pylin1991/nf-containers/optitype_1.3.5.sif \
    /tmp/optitype_1.3.5.def

# 確認版本
apptainer exec /data/pylin1991/nf-containers/optitype_1.3.5.sif \
    bash -c "OptiTypePipeline.py --help 2>&1 | head -3 \
             && razers3 --version 2>&1 | head -2 \
             && samtools --version | head -1"
```

---

## 資料庫下載與建置

### 目錄規劃

```
/data/pylin1991/GenomicReference/hg38/tertiary/    ← 原始下載
/scratch/pylin1991/GenomicReference_Cache/hg38/    ← nextflow 使用（ref_dir）
```

### hg38 reference FASTA

```bash
# 已在二級分析下載，直接用
ls /scratch/pylin1991/GenomicReference_Cache/hg38/Homo_sapiens_assembly38.fasta
```

### VEP cache 115（約 30 GB）

```bash
mkdir -p /data/pylin1991/GenomicReference/hg38/tertiary/vep_cache
cd /data/pylin1991/GenomicReference/hg38/tertiary/vep_cache
tmux new -s vep_download
wget -c https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz
tar xzf homo_sapiens_vep_115_GRCh38.tar.gz
rm homo_sapiens_vep_115_GRCh38.tar.gz
```

### dbNSFP 4.9c（含 P-KNN 整合）

```bash
# 手動下載：https://usf.box.com/shared/static/l8nik5s28i4zbup3b93hwz59dj2s94cp
mkdir -p /data/pylin1991/GenomicReference/hg38/tertiary/dbnsfp
cd /data/pylin1991/GenomicReference/hg38/tertiary/dbnsfp
unzip ~/Downloads/dbNSFP4.9c.zip

version=4.9c
zcat dbNSFP${version}_variant.chr1.gz | head -n1 > h
zgrep -h -v "^#chr" dbNSFP${version}_variant.chr*.gz \
    | sort -k1,1 -k2,2n -T /scratch/pylin1991/tmp \
    | cat h - \
    | bgzip -c > dbNSFP${version}_grch38.gz
tabix -s 1 -b 2 -e 2 dbNSFP${version}_grch38.gz
rm h dbNSFP4.9c_variant.chr*.gz

# 整合 P-KNN
conda activate genome
python3 /data/pylin1991/nf-containers/NGStertiary/1_0_0/scripts/build_dbnsfp_pknn.py \
    --dbnsfp   /scratch/pylin1991/GenomicReference_Cache/hg38/tertiary/dbnsfp/dbNSFP4.9c_grch38.gz \
    --pknn_dir /data/pylin1991/GenomicReference/hg38/tertiary/P_KNN_7 \
    --output   /scratch/pylin1991/GenomicReference_Cache/hg38/tertiary/dbnsfp/dbNSFP4.9c_with_pknn_grch38.gz
tabix -s 1 -b 2 -e 2 dbNSFP4.9c_with_pknn_grch38.gz
```

### dbNSFP 5.3a（`--academic_dbnsfp true` 用；含 P-KNN 整合）

只有開 `--academic_dbnsfp true` 才會用到。**兩份要並存**（4.9c 是預設的商用路徑，
5.3a 是學術路徑），不要互相覆蓋。設計理由與欄位差異見後面 §「dbNSFP 5.3a 與
--academic_dbnsfp」。

```bash
# 手動下載（與 4.9c 同一個發布頁；官方是 box.com 靜態連結，需登入頁面取得）
#   ⚠️ 實際用的連結請貼在這行下面，之後重建才追得回同一份
#   下載頁：https://sites.google.com/site/jpopgen/dbNSFP
cd /data/pylin1991/GenomicReference/hg38/tertiary/dbnsfp
unzip ~/Downloads/dbNSFP5.3a.zip

# 合併 per-chromosome 檔 → 單一 bgzip + tabix（流程與 4.9c 完全相同，只換 version）
#   若解壓後的檔名 pattern 不是 dbNSFP5.3a_variant.chr*.gz，改 version 變數即可
version=5.3a
zcat dbNSFP${version}_variant.chr1.gz | head -n1 > h
zgrep -h -v "^#chr" dbNSFP${version}_variant.chr*.gz \
    | sort -k1,1 -k2,2n -T /scratch/pylin1991/tmp \
    | cat h - \
    | bgzip -c > dbNSFP${version}_grch38.gz
tabix -s 1 -b 2 -e 2 dbNSFP${version}_grch38.gz
rm h dbNSFP5.3a_variant.chr*.gz

# 整合 P-KNN（在所有 dbNSFP 欄位後面附加一欄 PKNN_LLR）
conda activate genome
python3 /data/pylin1991/nf-containers/NGStertiary/1_0_0/scripts/build_dbnsfp_pknn.py \
    --dbnsfp   /scratch/pylin1991/GenomicReference_Cache/hg38/tertiary/dbnsfp/dbNSFP5.3a_grch38.gz \
    --pknn_dir /data/pylin1991/GenomicReference/hg38/tertiary/P_KNN_7 \
    --output   /scratch/pylin1991/GenomicReference_Cache/hg38/tertiary/dbnsfp/dbNSFP5.3a_with_pknn_grch38.gz
tabix -s 1 -b 2 -e 2 dbNSFP5.3a_with_pknn_grch38.gz
```

> ⚠️ **`sort` 不能省，而且要在 `cat h -` 之前。** `zgrep` 串接 `chr*.gz` 的順序是
> shell glob 的字典序（chr1, chr10, chr11 …），不是座標序。VEP 的 dbNSFP plugin 是
> **靠 tabix 隨機查詢**，index 建在未排序的檔案上會查不到（表現為分數欄大量變 `.`，
> 不會報錯）。4.9c 就是在這裡踩過 —— 合併後染色體順序不一致，得重新 `sort` →
> 重新 `bgzip` → 重新 `tabix`。`-T` 指到 scratch，不然 `/tmp` 會爆。
>
> ⚠️ **輸出寫在 scratch**（`/scratch/.../GenomicReference_Cache`）。scratch 不是永久
> 儲存，重建流程要能重跑；原始 zip 留在 `/data`。

#### 建置後驗證（不需跑 pipeline）

```bash
cd /scratch/pylin1991/GenomicReference_Cache/hg38/tertiary/dbnsfp
DB=dbNSFP5.3a_with_pknn_grch38.gz

# 1) tabix index 可用 + 座標真的排好了（隨機取一段查得到才算成功）
#    ⚠️ dbNSFP 的 #chr 欄不帶 chr 前綴（是 1/2/X），region 要照檔案裡的寫法
tabix -l $DB | head -3                           # 先看實際 contig 命名
CHR=$(tabix -l $DB | grep -x -m1 -e 17 -e chr17) # 兩種命名都接
tabix $DB ${CHR}:7676000-7676500 | wc -l         # TP53 區域，應 > 0
zcat $DB | cut -f1 | uniq | head -30             # 染色體應成塊且遞增，不是交錯

# 2) 確認 pipeline 要抓的欄位都存在（缺一個 VEP 就整欄變 "."，且不會報錯）
zcat $DB | head -1 | tr '\t' '\n' | grep -nE \
  '^(PKNN_LLR|REVEL_score|MutPred2_score|MutPred2_pred|VEST4_score|CADD_phred)$'
zcat $DB | head -1 | tr '\t' '\n' | grep -nE \
  '^gnomAD2\.1\.1_exomes_non_cancer_(AF|EAS_AF)$'   # ACMG 用的族群頻率（與 4.9c 對齊）
zcat $DB | head -1 | tr '\t' '\n' | grep -nE \
  '^gnomAD4\.1_joint_(AF|EAS_AF)$'                  # 參考欄，不進 ACMG 計分

# 3) PKNN_LLR 真的有值（最後一欄）
zcat $DB | awk -F'\t' 'NR>1 && $NF!="." {c++} END{print "PKNN_LLR 有值:", c}' | head -1

# 4) build_dbnsfp_pknn.py 的雙向統計（跑的時候印在 stderr，要留存）
#    重點看「P-KNN 載入但沒被用到」= 0；5.3a 實測零遺漏（P-KNN 原生就是 5.3 產生）
```

#### 測試（跑 pipeline，`--academic_dbnsfp true`）

```bash
nextflow -c .../nextflow_tertiary.config run .../main_tertiary.nf -profile local \
    --samplesheet samplesheet_nckuh.csv \
    --out_dir /scratch/pylin1991/tertiary_test_53a \
    --academic_dbnsfp true \
    -resume

# ⚠️ 換 dbNSFP 會讓 VEP_ANNOTATE 之後的所有步驟重算，-resume 只救得到 PREPARE_VCF。
#    想同時保留 4.9c 的結果，out_dir 要另開一個（不要蓋掉原本那份）。

# 檢查 1：banner 確認真的換了檔案（不是只有 flag 開著）
#   main_tertiary.nf 的 banner 印「實際使用」的那一份 + 模式，跑起來就在 console；
#   事後要查就 grep .nextflow.log（log.info 也會寫進去）
grep -E "dbNSFP( |\s)*(:|模式)" .nextflow.log | tail -4
#   期望：dbNSFP : ...dbNSFP5.3a_with_pknn_grch38.gz
#         dbNSFP 模式 : 5.3a（--academic_dbnsfp 啟用：...）
#   ※ preflight 現在也會檢查「實際使用的那一份」＋它的 .tbi，缺檔會直接 error 開不了跑

# 檢查 2：輸出的 DBNSFP_VERSION 欄
TSV=/scratch/pylin1991/tertiary_test_53a/<sample>/03_acmg/<sample>.snv_indel.acmg.tsv
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1{print $h["DBNSFP_VERSION"]}' $TSV \
    | sort | uniq -c            # 期望全部 5.3a

# 檢查 3：新增欄位的有值比例。**要在 missense SNV 母體內看**，
#         全表比例會被 synonymous/intron/UTR/indel 稀釋成 ~25%（dbNSFP 只收 nsSNV）
awk -F'\t' -v cols=REVEL_SCORE,CADD_PHRED,PKNN_LLR,MUTPRED2_SCORE,VEST4_SCORE,ALPHAMISSENSE_SCORE '
  NR==1 { for (i=1;i<=NF;i++) h[$i]=i
          n=split(cols,want,","); next }
  $h["CONSEQUENCE"] ~ /missense/ {
          tot++
          for (i=1;i<=n;i++) if ($h[want[i]] != "." && $h[want[i]] != "") k[want[i]]++ }
  END   { printf "missense SNV = %d\n", tot
          for (i=1;i<=n;i++) printf "  %-20s %6d  %5.1f%%\n", want[i], k[want[i]], 100*k[want[i]]/tot }
' $TSV
#   ⚠️ 欄名以輸出 TSV 的 header 為準（見附錄 A 的 81 欄總表）；
#      對不到會印 0%，那是欄名寫錯不是資料庫壞掉 —— 先用
#      `head -1 $TSV | tr '\t' '\n' | grep -n .` 對一下。
```

實測數字（NA12878 WES）見後面 §「驗證（NA12878 WES，`--academic_dbnsfp true`）」：
missense 母體內 REVEL 95.9%、CADD 99.5%、P-KNN 99.5%、MutPred2 98.9%、AlphaMissense 98.8%。
**用這組數字當回歸基準** —— 若重建後掉到明顯低於這些，先懷疑 `sort`／tabix，而不是資料庫本身。

### ClinVar

> ⚠️ NCBI 官方 ClinVar VCF contig 格式是 `1`, `2`（無 chr 前綴），必須 rename 否則 VEP annotation 全部為 `.`

```bash
cd /data/pylin1991/GenomicReference/hg38/tertiary/clinvar
wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_20260510.vcf.gz
wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_20260510.vcf.gz.tbi

for i in $(seq 1 22) X Y; do echo "$i chr$i"; done > /tmp/chr_rename.txt
echo "MT chrM" >> /tmp/chr_rename.txt

conda activate genome
bcftools annotate --rename-chrs /tmp/chr_rename.txt \
    clinvar_20260510.vcf.gz -Oz -o clinvar_20260510.chr.vcf.gz
tabix -p vcf clinvar_20260510.chr.vcf.gz
mv clinvar_20260510.chr.vcf.gz clinvar_20260510.vcf.gz
mv clinvar_20260510.chr.vcf.gz.tbi clinvar_20260510.vcf.gz.tbi

# 建 lookup table
wget "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2026-05.txt.gz"
python3 /data/pylin1991/nf-containers/NGStertiary/1_0_0/scripts/build_clinvar_lookup.py \
    --input  variant_summary_2026-05.txt.gz \
    --output /scratch/pylin1991/GenomicReference_Cache/hg38/tertiary/clinvar/clinvar_lookup.tsv.gz
rm variant_summary_2026-05.txt.gz
```

### gnomAD v4.1 genome（約 600 GB）

```bash
tmux new -s gnomad_download
mkdir -p /data/pylin1991/GenomicReference/hg38/tertiary/gnomad
cd /data/pylin1991/GenomicReference/hg38/tertiary/gnomad
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
    wget -c \
        "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz" \
        "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz.tbi"
done
```

### gnomAD v3.1 mito（CC0，取代 MITOMAP CC BY-NC）

```bash
mkdir -p /data/pylin1991/GenomicReference/hg38/tertiary/gnomad_mito
cd /data/pylin1991/GenomicReference/hg38/tertiary/gnomad_mito
wget -c "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"
wget -c "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz.tbi"

apptainer exec --bind /data,/scratch \
    /data/pylin1991/nf-containers/tertiary_python_1.0.0.sif \
    python3 /data/pylin1991/nf-containers/NGStertiary/1_0_0/scripts/build_gnomad_mito_lookup.py \
        --vcf    gnomad.genomes.v3.1.sites.chrM.vcf.bgz \
        --output gnomad_mito_lookup.tsv.gz
```

### LOFTEE 資料檔

```bash
mkdir -p /data/pylin1991/GenomicReference/hg38/tertiary/loftee
cd /data/pylin1991/GenomicReference/hg38/tertiary/loftee
wget -c https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz
wget -c https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai
wget -c https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi
wget -c https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz
wget -c https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw
gunzip loftee.sql.gz  # sql 需解壓，bw 和 fa.gz 不要解壓
curl -sL https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/115/LoFtool_scores.txt \
    -o LoFtool_scores.txt
```

### Pangolin annotation database

```bash
mkdir -p /data/pylin1991/GenomicReference/hg38/tertiary/pangolin
cd /data/pylin1991/GenomicReference/hg38/tertiary/pangolin
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz

apptainer exec --bind /scratch,/data \
    /data/pylin1991/nf-containers/pangolin_cu130_1.0.0.sif \
    create_db.py \
    --filter MANE_Select,MANE_Plus_Clinical,Ensembl_canonical \
    gencode.v47.annotation.gtf.gz
```

### STRchive + 查表

```bash
mkdir -p /data/pylin1991/GenomicReference/hg38/tertiary/strchive
cd /data/pylin1991/GenomicReference/hg38/tertiary/strchive
wget -c https://raw.githubusercontent.com/hdashnow/STRchive/main/data/STRchive-loci.json

conda activate genome
python3 /data/pylin1991/nf-containers/NGStertiary/1_0_0/scripts/build_str_lookup.py \
    --json         STRchive-loci.json \
    --output_varid str_lookup_varid.tsv.gz \
    --output_pos   str_lookup_pos.tsv.gz

# 複製到 ref cache
cp str_lookup_varid.tsv.gz str_lookup_pos.tsv.gz \
    /scratch/pylin1991/GenomicReference_Cache/hg38/tertiary/strchive/
```

### ClinGen（PVS1 + MOI）

```bash
mkdir -p /scratch/pylin1991/GenomicReference_Cache/hg38/tertiary/clingen
cd /scratch/pylin1991/GenomicReference_Cache/hg38/tertiary/clingen

wget https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
wget "https://search.clinicalgenome.org/kb/gene-validity/download" \
     -O clingen_gene_disease_validity.csv

conda activate genome
python3 /data/pylin1991/nf-containers/NGStertiary/1_0_0/scripts/build_gene_moi.py \
    --clingen_gene clingen_gene_disease_validity.csv \
    --omim_genemap /scratch/pylin1991/GenomicReference_Cache/hg38/tertiary/omim/genemap2.txt \
    --output gene_moi.tsv.gz
```

### AnnotSV Annotation Databases（約 2-3 GB）

```bash
tmux new -s annotsv_annotations
mkdir -p /data/pylin1991/GenomicReference/hg38/tertiary/annotsv_annotations
cd /data/pylin1991/GenomicReference/hg38/tertiary/annotsv_annotations

wget --tries=0 --wait=10 --retry-connrefused \
    https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.5.tar.gz
mkdir -p share/AnnotSV
tar xzf Annotations_Human_3.5.tar.gz -C share/AnnotSV/
rm Annotations_Human_3.5.tar.gz

# 驗證
apptainer exec \
    --bind /data/pylin1991/GenomicReference,/scratch/pylin1991 \
    /data/pylin1991/nf-containers/annotsv_3.5.10.sif \
    AnnotSV \
        -SVinputFile /scratch/pylin1991/Pipeline_test/NA12878_WES_PON/NA12878_WES/05_cnv_sv/NA12878_WES.delly.vcf.gz \
        -outputDir /scratch/pylin1991/annotsv_test \
        -annotationsDir share/AnnotSV \
        -genomeBuild GRCh38 \
        -SVinputInfo 1
ls /scratch/pylin1991/annotsv_test/
```


### StellarPGx repo 設定

```bash
cd /data/pylin1991/nf-containers/NGStertiary/1_0_0/
git clone https://github.com/SBIMB/StellarPGx stellarpgx_repo

# 確認 bgzip 格式（不同時間 clone 可能有格式問題）
conda activate genome
for f in stellarpgx_repo/resources/cyp2d6/res_hg38/*.vcf.gz; do
    if ! bgzip -t $f 2>/dev/null; then
        echo "修復：$f"
        mv $f ${f}.bak
        zcat ${f}.bak | bgzip -c > $f
        tabix -f -p vcf $f
    fi
done
```

> ⚠️ `stellarpgx_repo` 路徑在各環境的 config profile 裡各自設定（`local` / `dgm`），不要寫死在全域 params 裡，否則 Nextflow params 優先順序（global > profile）會讓 profile 的值被 null 蓋掉。

### PharmCAT positions VCF 建置

GATK HaplotypeCaller 用 `-L pharmcat_positions.vcf.gz` 只跑 PGx 位點，需要從 PharmCAT 容器取出：

```bash
mkdir -p /data/pylin1991/GenomicReference/hg38/tertiary/pharmcat

# 從 pharmcat_3.2.0.sif 取出（容器內建在 /pharmcat/pharmcat_positions.vcf）
apptainer exec /data/pylin1991/nf-containers/pharmcat_3.2.0.sif     cat /pharmcat/pharmcat_positions.vcf     > /data/pylin1991/GenomicReference/hg38/tertiary/pharmcat/pharmcat_positions.vcf

conda activate genome

# bgzip + tabix（必須 TBI format，GATK 不接受 CSI）
bgzip /data/pylin1991/GenomicReference/hg38/tertiary/pharmcat/pharmcat_positions.vcf
tabix -p vcf /data/pylin1991/GenomicReference/hg38/tertiary/pharmcat/pharmcat_positions.vcf.gz

# 用 GATK 建 TBI index（確保 GATK 可讀）
apptainer exec --bind /data,/scratch     /data/pylin1991/nf-containers/gatk_4.6.2.0.sif     gatk IndexFeatureFile         -I /data/pylin1991/GenomicReference/hg38/tertiary/pharmcat/pharmcat_positions.vcf.gz

# 確認位點數（應為 1207）
bcftools view -H     /data/pylin1991/GenomicReference/hg38/tertiary/pharmcat/pharmcat_positions.vcf.gz     | wc -l
```

> ⚠️ 不能直接用容器內的 `pharmcat_positions.vcf.bgz`（CSI index），GATK 只接受 TBI format。
> ⚠️ 需要同時有 GATK 產生的 `.tbi` 才能讓 GATK 正確處理所有 1207 個位點（否則只跑 chr1）。

---

## Pipeline 執行

### Sample sheet 格式

```csv
sample_id,pipeline_type,input_dir,seq_type,hpo
NA12878_WES,nckuh,/scratch/pylin1991/Pipeline_test/NA12878_WES_PON/NA12878_WES/NA12878_WES,WES,
NA12878_WGS,nckuh,/scratch/pylin1991/Pipeline_test/NA12878_WGS/NA12878_WGS/NA12878_WGS,WGS,
VAL-10,dragen,/scratch/pylin1991/Pipeline_test/DRAGEN/VAL-10,WGS,HP:0001250
```

```csv
sample_id,pipeline_type,input_dir,seq_type,hpo
26T00001,dragen,/data/WGS/26T00001,WGS,
26T00089,dragen,/data/WGS/26T00089,WGS,
26T00076,dragen,/data/WGS/26T00076,WGS,
```


    

### 執行指令

```bash
conda activate nextflow
cd /data/pylin1991/nf-containers/NGStertiary/1_0_0/

# NCKUH（WES + WGS 混合，PGx 預設開啟）
nextflow -c nextflow_tertiary.config run main_tertiary.nf \
    -profile local --pipeline_type nckuh \
    --samplesheet /scratch/pylin1991/tertiary_test/samplesheet_nckuh.csv \
    --out_dir /scratch/pylin1991/tertiary_test \
    -resume

# DRAGEN
nextflow -c nextflow_tertiary.config run main_tertiary.nf \
    -profile local --pipeline_type dragen \
    --samplesheet /scratch/pylin1991/tertiary_test/samplesheet_nckuh.csv \
    --out_dir /scratch/pylin1991/tertiary_test \
    --run_phasing true \
    -resume

# 關閉 PGx（速度測試用）
nextflow -c nextflow_tertiary.config run main_tertiary.nf \
    -profile local --pipeline_type nckuh \
    --samplesheet /scratch/pylin1991/tertiary_test/samplesheet_nckuh.csv \
    --out_dir /scratch/pylin1991/tertiary_test \
    --run_pgx false \
    -resume

# DGM Server
nextflow -c nextflow_tertiary.config run main_tertiary.nf \
    -profile dgm --pipeline_type nckuh \
    --samplesheet /home/pipeline/samplesheet.csv \
    --out_dir /home/pipeline/tertiary_output \
    -resume
```

---

## 結果驗證

### 驗證指令（NA12878）

```bash
BASE=/scratch/pylin1991/tertiary_test
SAMPLE=NA12878_WES   # 或 NA12878_WGS

# SNV/Indel 行數
wc -l $BASE/$SAMPLE/03_acmg/${SAMPLE}.snv_indel.acmg.tsv

# ACMG 分類分布
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="ACMG_CLASS") col=i}
            NR>1{print $col}' \
    $BASE/$SAMPLE/03_acmg/${SAMPLE}.snv_indel.acmg.tsv \
    | sort | uniq -c | sort -rn

# MITO
wc -l $BASE/$SAMPLE/04_mito/${SAMPLE}.mito.tsv

# STR
wc -l $BASE/$SAMPLE/05_str/${SAMPLE}.str.tsv

# CNV/SV
ls -lh $BASE/$SAMPLE/06_cnv_sv/

# PGx
ls -lh $BASE/$SAMPLE/07_pgx/
cat $BASE/$SAMPLE/07_pgx/${SAMPLE}.optitype.tsv
grep "CYP2D6\|CYP2C9\|HLA\|GENE" \
    $BASE/$SAMPLE/07_pgx/${SAMPLE}.pgx.tsv \
    | cut -f3-6 | sort -u
```

### 預期結果（NA12878）

| Module | NA12878_WES | NA12878_WGS |
|---|---|---|
| SNV/Indel | ~37,199 行 | ~5,729,808 行 |
| ACMG LP | ~41（低 VAF artifact）| ~12 |
| MITO | ~33 行 | ~21 行 |
| STR | ~6 行 | ~19 行 |
| CNV | ~109K | ~4.9M |
| SV | ~12M | ~820M / 957K 行 |
| PGx pgx.tsv 行數 | ~136 行 | ~186 行 |
| PGx Called 基因數 | 17/23 | 22/23 |
| PGx CYP2D6 | VCF-only（WES，準確度低）| `*1/*5`, activity=1.0 ✅ |
| PGx HLA-A | N/A（WES 無 BAM）| `*01:01/*01:01` ✅ |
| PGx HLA-B | N/A（WES 無 BAM）| `*08:01/*08:01`（ground truth `*07:02/*40:02`）|
| PGx MT-RNR1 | Reference（LOW risk）✅ | Reference（LOW risk）✅ |
| PGx 臨床警示 | N/A | `*57:01 negative`, `*58:01 negative` ✅ |

> HLA-B heterozygous allele 辨認是 OptiType 的已知限制，在 reads 不夠多時容易 call 成 homozygous。臨床最重要的 `*57:01`（abacavir）和 `*58:01`（allopurinol）negative 結果正確，臨床安全性無虞。


---

## 傳送至 DGM Server

### 傳送容器
```bash
# vep_115.sif（約 5GB）
# pangolin_cu121_1.0.0.sif（DGM 用這顆，cu121 —— 與 DGX-2 同一顆）
#   DGM 那張卡兩顆都跑得動（見容器建立章節的機器表）；統一用 cu121 是為了
#   評鑑能聲明「所有 production 使用同一顆容器」，見容器建立章節的說明。
#   ※ pangolin_cu130_1.0.0.sif（cu130）只有開發機需要，不必傳到 DGM
# tertiary_python_1.0.0.sif
scp /data/pylin1991/nf-containers/*.sif \
    n101569@192.168.84.91:/home/pipeline/nextflow_containers/

rsync -avz --progress \
    /data/pylin1991/nf-containers/*.sif \
    n101569@192.168.84.91:/home/pipeline/nextflow_containers/

```

### 傳送 Reference 資料庫（rsync，支援斷點續傳）
```bash
# 不傳 gnomad（約 600GB，VEP cache 內建版本已足夠 Phase 1）
rsync -avz --progress \
    --exclude='gnomad/' \
    /data/pylin1991/GenomicReference/hg38/ \
    n101569@192.168.84.91:/home/pipeline/reference/hg38/
```

### 傳送 Pipeline 程式碼
```bash
rsync -avz --progress \
    /data/pylin1991/nf-containers/NGStertiary/1_0_0/ \
    n101569@192.168.84.91:/home/pipeline/tertiary_code/

scp -r /data/pylin1991/nf-containers/NGStertiary/1_0_0/ n101569@192.168.84.91:/home/pipeline/tertiary_code/
```

### DGM 執行
```bash
ssh n101569@192.168.84.91

source /home/pipeline/pipeline_code/DGM_NGS2ndAnalysis.sh

nextflow -c /home/pipeline/tertiary_code/nextflow_tertiary.config \
    run /home/pipeline/tertiary_code/main_tertiary.nf \
    -profile dgm \
    --samplesheet /home/pipeline/samplesheet_nckuh.csv \
    --out_dir /home/pipeline/tertiary_output \
    -resume
```

> ⚠️ v3.1 起改為 sample sheet 批次輸入，舊的 `--sample_id / --input_dir / --seq_type`
> 三個參數已不再使用。

---

## 傳送至 DGX-2

DGX-2 與二級分析共用同一個 reference 與容器目錄（`/datalake_Intermediate/pipeline/`），
但**三級程式碼放獨立的 `tertiary_code/`**，不與二級的 `pipeline_code/` 混用（兩者都有
`modules/` 和 `scripts/`，放同一層會互相覆蓋）。

| 項目 | 位置 |
|------|------|
| 帳號 | `n101569@10.11.33.75` |
| Reference | `/datalake_Intermediate/pipeline/reference/hg38`（**與二級共用**）|
| 容器 | `/datalake_Intermediate/pipeline/nextflow_containers`（**與二級共用**）|
| 三級程式碼 | `/datalake_Intermediate/pipeline/tertiary_code` |
| 輸出 | `/datalake_Intermediate/pipeline/nextflow_output` |
| GPU | V100 × 6，與二級共用 → **必須用 GPU lock**（`-profile dgx` 已預設開啟）|

### 建立資料夾（僅三級需要新增的部分）
```bash
ssh n101569@10.11.33.75
mkdir -p /datalake_Intermediate/pipeline/tertiary_code
# reference / nextflow_containers / nextflow_output 由二級部署時已建立
```

### 傳送容器（三級專用的幾個；其餘與二級共用）
```bash
rsync -avz --progress \
    /data/pylin1991/nf-containers/{vep_115,pangolin_cu121_1.0.0,tertiary_python_1.0.0,annotsv_3.5.10,pharmcat_3.2.0,stellarpgx_graphtyper2.5.1,optitype_1.3.5,samtools_1.23.1}.sif \
    n101569@10.11.33.75:/datalake_Intermediate/pipeline/nextflow_containers/
```

### 傳送 Reference（三級資料庫）
```bash
# 只傳 tertiary/ 子目錄即可（hg38 主參考二級已傳過）
# gnomad/ 約 600GB，VEP cache 內建版本已足夠 → 排除
rsync -avz --progress --exclude='gnomad/' \
    /data/pylin1991/GenomicReference/hg38/tertiary/ \
    n101569@10.11.33.75:/datalake_Intermediate/pipeline/reference/hg38/tertiary/
```

### 傳送 Pipeline 程式碼
```bash
rsync -avz --progress \
    /data/pylin1991/nf-containers/NGStertiary/1_0_0/ \
    n101569@10.11.33.75:/datalake_Intermediate/pipeline/tertiary_code/
```

### DGX-2 執行
```bash
ssh n101569@10.11.33.75

nextflow -c /datalake_Intermediate/pipeline/tertiary_code/nextflow_tertiary.config \
    run /datalake_Intermediate/pipeline/tertiary_code/main_tertiary.nf \
    -profile dgx \
    --samplesheet /datalake_Intermediate/pipeline/samplesheet_nckuh.csv \
    --out_dir /datalake_Intermediate/pipeline/nextflow_output \
    -resume
```

### 部署後檢查
```bash
# 1) 路徑推導是否正確（不需真樣本）
nextflow -c .../nextflow_tertiary.config run .../main_tertiary.nf -profile dgx \
    --samplesheet /dev/null --out_dir /tmp/x -preview

# 2) Pangolin 能否吃 V100（compute 7.0 / sm_70）
SIF=/datalake_Intermediate/pipeline/nextflow_containers/pangolin_cu121_1.0.0.sif
apptainer exec --nv $SIF \
    python3 -c "import torch; print('cuda', torch.cuda.is_available(), torch.cuda.get_device_name(0))"
#   期望：cuda True Tesla V100-SXM3-32GB
#
#   ⚠️ is_available()=True 還不夠！arch 不符時它照樣 True，要跑到 kernel 才爆。
#      一定要再跑一次真的 forward：
apptainer exec --nv $SIF python3 -c \
    "import torch; x=torch.randn(1,4,64,device='cuda'); \
     c=torch.nn.Conv1d(4,8,3).cuda(); print('conv OK', c(x).shape)"
#
#   ⚠️ 兩種錯誤都是「容器的 torch wheel 不對」，不是 V100 不行、也不是驅動該升：
#      - "driver ... too old (found version 12020)"        → wheel 是 CUDA 13.x
#      - "no kernel image is available ... on the device"  → wheel 缺這張卡的 arch
#      先查容器的 CUDA 版本與 arch flags：
apptainer exec $SIF python3 -c \
    "import torch; print(torch.__version__, torch.version.cuda); \
     print(torch._C._cuda_getArchFlags())"
#      arch flags 缺 sm_70 → 這顆不是 cu121 版（很可能誤傳了開發機的 cu130 容器）。
#      production 只能用 pangolin_cu121_1.0.0.sif；第一行 TARGET 可以確認：
apptainer exec $SIF head -1 /opt/build_versions.txt      # 期望 TARGET=cu121
#      見「容器建立 → Pangolin → 為什麼是兩顆容器」。
#      重建期間可先用 --use_gpu_pangolin false 走 CPU（結果相同，只影響速度）。

# 3) GPU lock 腳本存在（與二級共用）
ls -l /datalake_Intermediate/pipeline/pipeline_code/gpu_{lock,unlock}.sh
```

> **GPU lock**：DGX-2 是共用機器，二級（Parabricks）與三級（Pangolin）可能同時執行。
> `-profile dgx` 會啟用 `use_gpu_lock`，每個 Pangolin task 由 `gpu_lock.sh` 搶一張空閒卡、
> `trap EXIT` 歸還；`maxForks = 6` 對應六張 V100。**卡號不要寫死在 config**，否則多個 task
> 會擠同一張。

---

## PGx Module 建置記錄

### 工具與容器

| 工具 | 版本 | License | 容器 | 用途 |
|---|---|---|---|---|
| PharmCAT | 3.2.0 | MPL 2.0 ✅ | `pharmcat_3.2.0.sif` | star allele calling + CPIC/DPWG |
| StellarPGx | 1.2.8（graphtyper 2.5.1）| Open source ✅ | `stellarpgx_graphtyper2.5.1.sif` | CYP2D6 BAM-based（WGS only）|
| OptiType | 1.3.5 | BSD-3-Clause ✅ | `optitype_1.3.5.sif`（自建）| HLA-A/B typing（WGS only）|
| SAMtools | 1.23.1 | MIT ✅ | `samtools_1.23.1.sif` | HLA reads 擷取 |

> ❌ 禁用：BCyrius（PolyForm Strict）、Aldy（non-commercial）

### PGx workflow

```
WGS + WES 所有樣本（有 BAM）：
  BAM → PGX_GVCF（GATK HaplotypeCaller，只跑 1207 個 PGx 位點）
      → gVCF（含 0/0 reference call）→ PGX_PHARMCAT
  BAM → PGX_MTRN1（bcftools mpileup，chrM:827/1494/1555，< 1 秒）
      → mtrn1.vcf.gz → PGX_PARSE（MT-RNR1 coverage 確認）

WGS only：
  BAM → PGX_STELLARPGX（StellarPGx）→ CYP2D6 diplotype（outside call）
  BAM → PGX_HLA_EXTRACT（samtools）  → chr6 region + 16 alt contigs → fastq
      → PGX_OPTITYPE（razers3 + OptiType）→ HLA-A/B allele（outside call）

PharmCAT：
  gVCF + outside_calls → pharmcat_vcf_preprocessor → pharmcat.jar → report.json

PGX_PARSE：
  report.json + mito.tsv + mtrn1.vcf.gz → parse_pgx_report.py → pgx.tsv（16 欄）
```

**gVCF 的重要性：**
ensemble/DRAGEN VCF 只有 variant calls，PharmCAT 無法區分 reference（0/0）vs no-coverage（missing）。GATK gVCF 讓 PharmCAT 幾乎所有基因都能 call 出正確 diplotype：

| | gVCF 前 | gVCF 後 |
|---|---|---|
| WGS Called 基因數 | ~8 個 | **22 個** |
| WGS Unknown | ~14 個 | **1 個**（MT-RNR1，已由 mito pipeline 補充）|
| WES Called 基因數 | ~5 個 | **17 個** |

### HLA reads 擷取策略

hg38 的 HLA reads 分散在三處：
1. `chr6:29940260-33086201`（HLA-A/B/C/DRB1 核心區域）
2. 16 個 chr6 alt contigs（`chr6_GL000250v2_alt` 等，HLA alt haplotype）

**不取 unmapped reads（最終決策）：**
- NCKUH BWA：unmapped 只有 11.4 萬（chr6 region 181 萬，貢獻 < 6%）
- DRAGEN：unmapped 高達 3000 萬（DRAGEN 保留所有無法比對的 reads），若取 unmapped 導致 razers3 跑 5+ 小時
- 亞洲 HLA allele（B*46:01, B*58:01）在 hg38 primary assembly 正常比對到 chr6，不需要 unmapped
- 統一只取 chr6 region + alt contigs，兩個 pipeline 行為一致

**為何 single-end**：從 BAM 重建的 reads 配對關係不完整。OptiType paired-end 模式預設 `unpaired_weight=0`，強行拆 R1/R2 會丟棄大量 singleton（只剩 ~44 reads）。Single-end 模式所有 reads 都計分，實測 reads 數量更多（~6,475 reads）。

**razers3 -i 90**：預設 97% identity 過濾太嚴（~78 reads），`-i 90` 讓更多 HLA reads 通過。

### 與其他 pipeline 比較

| | 本 pipeline | 同學（BCyrius）|
|---|---|---|
| CYP2D6 方法 | Graphtyper2 graph re-genotyping | sequence matching |
| Gene deletion *5 | ✅ depth-based | ✅ |
| License | Open source ✅ | **PolyForm Strict ❌** |
| 速度 | 慢（Graphtyper2 重）| 快 |

| | 本 pipeline（OptiType）| 同學（OptiType）|
|---|---|---|
| HLA region | chr6 精確 region + alt contigs + unmapped | chr6 region only |
| razers3 過濾 | ✅ -i 90 | ❌ 跳過 |
| reads 模式 | single-end（從 BAM 重建）| paired-end |

---

## 踩雷記錄

### Apptainer

- mksquashfs 4.7.5 bug → `--disable-cache` + `APPTAINER_SQUASH_OPTIONS="-processors 1"` 兩個都要加
- `conda activate nextflow` 後跑 `apptainer build` 會失敗，要用 `conda activate base`

### Nextflow config params 優先順序

- **global params > profile params**：不要在全域 params 宣告 `stellarpgx_repo = null`，否則會蓋掉 profile 的值
- 正確做法：只在 profile 裡設定，全域不宣告

### PharmCAT 3.2.0

- `pharmcat_pipeline` 沒有 `-po` 參數 → 必須拆成 `pharmcat_vcf_preprocessor` + `pharmcat.jar -po`
- `pharmcat_vcf_preprocessor` 在 3.x 是編譯好的執行檔，不是 Python script
- NCKUH VCF 有雙 sample column（`_DV` + `_HC`）→ `bcftools query -l | head -1` 取第一個
- `activityScore` 在 JSON 裡存在但可能是 `null`（大多數基因本來就無 activity score）
- RECOMMENDATION 有 HTML entities → `strip_html()` 需加 entity decode

### StellarPGx

- GitHub repo 的 `.sif` 是 Git LFS 指標 → 必須從 Docker Hub pull
- `stellarpgx_repo` 的 `.vcf.gz` 可能是普通 gzip 不是 bgzip → 需重新用 bgzip 壓縮
- 不同環境 clone 的 repo 可能格式不同 → 每次 clone 後用 `bgzip -t` 確認

### OptiType

- `fred2/optitype:latest` → Python 2.7/3.5，`parse_optitype.py` 需相容：
  - f-string 改 `.format()`
  - `list[dict]` type hint 改成無 type annotation
- razers3 舊版不支援 `.fastq` 輸出（只有 `.bam`, `.fasta`, `.sam`）
- razers3 paired-end 模式記憶體暴增 → `std::bad_alloc`（exit 134）
- OptiType 1.5.0 的 Pyomo 6.10 有 constraint infeasible bug（reads 少時）
- razers3 `--no-home` + `--home /tmp` 才能解決 `/home/{user}` mkdir 問題

### Nextflow channel 設計

- WES 沒有 BAM → 不能用 `join/remainder` 策略
- 正確做法：按 `seq_type` 分成 `pgx_wgs_ch` 和 `pgx_wes_ch`，各自進 `PGX_ANNOTATE`
- `pipeline_type` 變數要從 `samples` 正確取得（不能用未定義的全域變數）

### GATK HaplotypeCaller（PharmCAT gVCF）

- pharmcat_positions.vcf.bgz（從容器複製）的 index 是 CSI format → GATK 不接受
- 必須重新 `cat /pharmcat/pharmcat_positions.vcf | bgzip > pharmcat_positions.vcf.gz && tabix` + `gatk IndexFeatureFile` 產生 TBI format
- `-L pharmcat_positions.vcf.gz`（VCF 當 interval）只跑了 chr1（24,280 bp），其他染色體沒跑到 → 必須同時用 `gatk IndexFeatureFile` 建 TBI 才能讓 GATK 正確讀取所有 1207 個位點
- GATK `--output-mode EMIT_ALL_ACTIVE_SITES` 輸出的 0/0 位點 QUAL 欄位是 `inf` → PharmCAT VCF parser 不接受 → 用 `awk` 換成 `99`

### Transcript Picking（parse_vep_csq.py）

- **問題**：舊邏輯把所有基因的 transcript 混在一起，鄰近基因的 MANE Select（如 TRAPPC4 downstream）蓋掉本基因的 HIGH impact transcript（如 SLC37A4 stop_gained/Pathogenic）
- **根本原因**：SLC37A4 的 MANE Select（ENST00000642844）在 alt contig NW_009646203.1 上，主染色體 VEP annotation 不會出現，所以 SLC37A4 所有主染色體 transcript 都沒有 MANE_SELECT 標記
- **MANE summary 確認**：`MANE.GRCh38.v1.5.summary.txt.gz` 中 SLC37A4 的 MANE Select 是 NM_001164277.2 / ENST00000642844，對應的 contig 是 NW_009646203.1
- **修法（v3.2）**：`pick_transcripts_for_output()` 按基因分組，每個基因選代表 transcript，再從所有代表中保留所有 MANE 和比 MANE 更嚴重的 non-MANE
- **MANE_ALL 欄位移除**：新邏輯每個 transcript 各自一行，不需要 JSON 彙總
- **輸出欄位數**：61 欄 → 60 欄（移除 MANE_ALL）→ 61 欄（v3.5 新增 STRAND_BIAS，欄 24）

### Transcript Picking（parse_mito_vcf.py）

- **問題**：`pick_transcript` 只看 PICK=1，MT-CO1 upstream（protein_coding canonical）蓋掉 MT-TL1 non_coding_transcript_exon（Mt_tRNA）
- **測試 case**：chrM:3243 A>G（m.3243A>G，MELAS 最常見致病位點）→ 錯誤輸出 MT-CO1，正確應為 MT-TL1
- **修法**：優先選實際落在 gene body 內的 consequence，upstream/downstream/intergenic 降到最後（rank 99/100）
- **驗證**：VAL-33-WGS mito VCF 跑 parse_mito_vcf.py → chrM:3243 正確輸出 MT-TL1 ENST00000386347.1:n.14A>G

### MT-RNR1 mpileup

- mito.tsv 只有 variant calls，reference 位點不出現 → `parse_mito_tsv` 找不到任何位點 → 回傳空 → pgx.tsv 沒有 MT-RNR1
- 解法：`PGX_MTRN1`（bcftools mpileup，只跑 chrM:827/1494/1555）確認 coverage，有 coverage 才補 Reference
- bcftools mpileup `-r` region 必須按 coordinate 升序（827 < 1494 < 1555），否則 `hts_idx_push` 報錯
- bcftools mpileup 輸出後需加 `bcftools sort` 再 tabix，避免排序問題

#### ⚠️ coverage 閘門失效（2026-08 修正）

寫評鑑大補帖、逐行核對 `parse_pgx_report.py` 時抓到：上面那句「有 coverage 才補 Reference」
**當時並沒有真的生效**。

```python
if covered_positions or no_mito is False:   # ← no_mito 在這裡必然是 False
```

`no_mito == True` 的情況函式在前面就 `return` 了，所以 `no_mito is False` 恆為真，
`or` 短路掉整個 coverage 判斷 —— 只要 `mito.tsv` 存在且沒 call 到 MT-RNR1 致病變異，
**即使三個位點深度都是 0，也照樣輸出 `Reference / MTRN1_RISK=LOW`**，
`NOTES` 還會寫 "Coverage confirmed via mito pipeline"。
臨床上等於對一個沒測到的樣本宣告「可正常使用 aminoglycoside」→ 不可逆聽力損傷風險。

修正三處（`scripts/parse_pgx_report.py`）：

1. 拿掉 `or no_mito is False`，coverage 閘門真的生效；三點都 `DP < 10` → 不輸出列 → Unknown
2. coverage 只認 `MTRN1_PATHOGENIC` 的三個位點（原本走訪 VCF 全部位點，混進其他 chrM 位點會誤判）
3. `RECOMMENDATION` 只列**實際評估過**的變異；部分覆蓋時明列哪些位點未評估
4. `DP >= 10` 提為模組常數 `MTRN1_MIN_DP`

`MTRN1_RISK` 值域（`HIGH`/`LOW`/Unknown）不變，下游 GUI 不受影響。
7 個情境實測（全覆蓋 / 全 0 / 部分覆蓋 / 無 VCF / 只有其他位點 / 帶 m.1555A>G / 無 mito.tsv）皆符合預期。
詳見 `docs/評鑑大補帖_三級分析.md` §10.12。

### ClinVar

- NCBI 下載的 VCF contig 格式是 `1`, `2`（無 chr 前綴）→ VEP annotation 全部為 `.`
- 必須用 `bcftools annotate --rename-chrs` 轉換後才能使用

### LOFTEE gerp

- LOFTEE 的 `gerp_dist.pl` 有 bug，`loftee_path` 直接與檔名拼接（無 `/`）
- 解法：用 `--bind` 把 gerp bw 直接掛載至 `/opt/vep/Plugins/gerp_...bw`

---

## v3.5 更新記錄

### Sub-workflow 重構（main_tertiary.nf）

把 CNV/SV、STR、SNV 尾段的「flat channel + 裸 process」鏈包成 sub-workflow，與二級一致：

- `ANNOTATE_CNV_SV_NCKUH`（gCNV[WES] / CNVkit→BED[WGS] + Delly SV）、`ANNOTATE_CNV_SV_DRAGEN`
  （filter → AnnotSV）→ 皆定義在 `cnv_sv_annotation.nf`
- `ANNOTATE_STR_NCKUH`（PREPARE→PARSE）、`ANNOTATE_STR_DRAGEN`（單 process，對稱包起來）→ `str_annotation.nf`
- `ANNOTATE_SNV`（VEP→CSQ→ACMG，NCKUH/DRAGEN 共用）→ 新檔 `annotate_snv.nf`（跨三個 module，故獨立成檔）
- NCKUH CNV 的 WES/WGS `if (seq_types.contains(...))` 條件移除：改成無條件建 channel，空 channel = 0 task
  （語意相同，DAG 更完整）。`MITO_ANNOTATE`、`PLOIDY_REPORT_DRAGEN` 單 process 維持裸呼叫。

### STRAND_BIAS 欄位（parse_vep_csq.py，欄 24）

- 依 INFO/FS（FisherStrand）+ INFO/SOR 判股偏；germline 只標記不硬刪。
- `PASS` / `WARN(FS=..,SOR=..)`（GATK 門檻 SNV FS>60/SOR>3.0、indel FS>200/SOR>10.0）；
  `.` = 無 FS/SOR（DeepVariant-only 位點）→ 人工複核。DRAGEN/HC 位點有 FS/SOR。

### DRAGEN combined-record「AD 消失」修復（combine_phased.py）

- **症狀**：DRAGEN 走 `COMBINE_DRAGEN`（combine_phased.py，用原生 PS）合併 cis compound 後，combined MNV
  只剩 `GT:PS`，`AD/DP/VAF` 全掉成 `.`（VAL-58 chr17:80260571、VAL-10 145k+ 筆）。
- **修法**：combined record 改為**繼承 anchor（cluster 內最寬的 biallelic record）的完整 FORMAT**，只覆寫
  `GT`/`PS`；四種情況（重建出 2-ALT、無 biallelic anchor、mixed ploidy、chrM haploid）直接 passthrough。
- **驗證（2026-07，結案）**：VAL-10 重跑 → `03_acmg` 的 `AD_DV` 5,940,465 / 5,940,563 有值（僅 98 筆
  來源資料本身無 AD）。二級 `scripts/combine_phased.py` 與本 repo 逐位元組相同，md5 需一致。

### DRAGEN ploidy QC（PLOIDY_REPORT_DRAGEN + parse_dragen_ploidy.py）

- 讀 DRAGEN 原生 `{sample}.ploidy.vcf.gz` → `00_prepare/{sample}.ploidy_qc.txt`
  （性別核型 + 每 contig NDC + aneuploidy 警示）；找不到 ploidy.vcf 則 warn 後跳過。
- NDC 語意與二級 mosdepth `ploidy_check.py` 統一（正規化到估計核型的期望，~1.0 = 正常）。
- 驗證：VAL-10 → estimated `XX` / `sex_check: OK`。

### DRAGEN PGx 交叉註記（compare_dragen_pgx.py + PGX_DRAGEN_CONCORDANCE，僅 DRAGEN）

- 目的：把 DRAGEN 原生 PGx 判讀（`other/{sample}/germline_seq/{sample}.targeted.json` 的
  `locusAnnotations` + `cyp2d6`/`cyp2b6`）交叉比對到我們的 `pgx.tsv`，寫進 `NOTES`（**欄位不變**）：
  `DRAGEN 一致 / 不一致 / 未比對: <DRAGEN 原始 genotype>`。
- 正規化：reference-like（`Reference`/`wildtype`/`ref`/`*1`/`B(reference)`/`B(wildtype)`）→ `REF`；
  star allele 取排序集合；DRAGEN `;` 模糊多重解 → 我們的 diplotype 落在任一候選即「一致」；
  命名系統不同（star vs HGVS/rs）→ `未比對`（不妄下判定）。DRAGEN 沒有／我們沒有的基因不動。
- 位置：`PGX_DRAGEN_CONCORDANCE` 收在 `PGX_ANNOTATE` sub-workflow 內（`dragen_targeted_ch` 第 4 個
  take，非 DRAGEN 傳空 channel → 0 task）；為 `PGX_PARSE` 嚴格下游，發布的 pgx.tsv 取代 base 版。
- 驗證（VAL-10）：一致 10（含 UGT1A1 模糊候選、G6PD `B(reference)`↔`B(wildtype)`、CACNA1S/RYR1/
  MT-RNR1 reference）、不一致 2（CYP2D6 結構型排列、DPYD 我們 Indeterminate vs DRAGEN 解出 `*6`）、
  HLA-A/B 因 DRAGEN 無此基因不動。


---

## v3.6 更新記錄

### ClinVar 版本同步（config）

`params.clinvar`（VEP `--custom`，供 CLNSIG/CLNREVSTAT/CLNDN/CLNSIGCONF）與
`params.clinvar_lookup_tsv`（`parse_vep_csq.py` 查 Variation ID）**必須是同一個 ClinVar
release**。先前只重建了 lookup（2026-07），VCF 還停在 20260510 → 同一變異的致病性與 ID 會來自
不同版本。已改指 `clinvar_20260720.vcf.gz`，並在 config 註記此規則與 chr 改名＋tabix 的前提。

### ClinGen Evidence Repository 對照（build_clingen_erepo_lookup.py）

ERepo = 各 VCEP 專家小組的變異判讀，含**實際套用的 ACMG criteria**。以 ClinVar Variation ID
對照進我們的表，新增 `CLINGEN_VCEP_CLASS/_CRITERIA/_PANEL` 與 `CLINGEN_AGREEMENT`
（`AGREE` / `DIFFER_TIER` / `DIFFER` / `.`）。

- **只作對照、不進計分**：ClinGen SVI 2018 建議不要用 PP5/BP6（拿他人判讀當證據 = 循環論證），
  本 pipeline 亦未實作 PP5/BP6。
- 下載：`https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download`（TSV，20 欄）
- 實測 13,039 筆 → 輸出 12,852 個變異（98.6%）。
- **踩雷 1｜分隔符判斷**：原本從前 4KB 判斷 tab/comma，但 `HGVS Expressions` 欄一列就有 20+ 個
  逗號分隔的 HGVS → 被誤判成 CSV，整個 header 塌成一欄。改成**只看 header 那一行**。
- **踩雷 2｜Retracted**：ERepo 有已撤回的判讀（實測 8 筆），拿來當對照基準會誤導 → 過濾掉。
- **踩雷 3｜5% 沒有 ClinVar ID**：只有 Allele Registry ID。改用 HGVS 裡的 **GRCh38** accession
  （同列並存 NCBI36/GRCh37/GRCh38，必須指名 GRCh38）組 `chr:pos:ref:alt` 當備援 key，
  救回 511 筆；del/dup 的 g. 寫法無明確 REF/ALT，無法救（169 筆）。

### dbNSFP 5.3a 與 --academic_dbnsfp

- `--academic_dbnsfp true` → VEP 改用 5.3a，並多抓 REVEL / MutPred2 / VEST4 / CADD_phred。
  這些多為「學術免費、商業需授權」（CADD 尤其明確），故不放預設路徑。
- **不用 tabix 事後查表**：dbNSFP 分數欄是 `;` 分隔的多轉錄本值，自己查表等於重寫 VEP plugin
  的轉錄本配對邏輯，風險高；**也不跑兩輪 VEP**（VEP 是最慢的一步）。改為單次 VEP 切換檔案。
- **ACMG 基準不動**：兩版都用 gnomAD **2.1.1** exomes（4.9c 為整體；5.3a 只有子集 → 取
  `non_cancer`，~118k 最接近整體 ~125k）。gnomAD 4.1 另開 `GNOMAD41_JOINT_*` 參考欄，不計分。
  影響面實測很窄：`GNOMAD_E_AF_DBNSFP` 根本沒進 ACMG（純顯示），`GNOMAD_E_EAS_AF_DBNSFP`
  只在 AR/XL 的 PM2 作為 `min_eas_af()` 四個來源之一。
- 欄位取聯集且**附加在最後**，GUI schema 固定；`DBNSFP_VERSION` 記錄實際版本。

### P-KNN 合併：CSV 引號 bug（build_dbnsfp_pknn.py）

- **症狀**：P-KNN 有值、dbNSFP 也有完全相同的 pos/ref/alt，合併後卻是 `.`（chr1 實測 314 筆）。
- **原因**：部分 P-KNN 行有「引號包住、內含逗號」的欄位（如 `"Pathogenic,_no_conflicts"`），
  `line.split(",")` 會在引號內切開 → 欄位整體位移 → 第 21 欄（LLR）拿到 `_no_conflicts"` →
  `float()` 失敗 → **靜默跳過**。改用 `csv.reader`（處理引號）並放寬 `field_size_limit`。
- **雙向統計**：每條染色體印出「dbNSFP 有/無 LLR」＋「P-KNN 載入/被用到/**沒被用到**」。
  沒被用到 > 0 才是 key 對不上的警訊。無 LLR 的列另依 aaref/aaalt 分成 nonsense / stoploss /
  non-coding·splice / synonymous / **missense**——P-KNN 是 missense-only，只有 missense 對不到
  才是真問題。
- **實測結論**：5.3a 零遺漏；4.9c 較差是因為 P-KNN 原生就是 5.3 產生（版本落差），非 bug。
  剩下對不到的 missense 約 1.8%，是 P-KNN 只用 **MANE Select** 產生所致（設計選擇）。

### PVS1 改為 ClinGen SVI 決策樹（Abou Tayoun 2018）

舊版：LOFTEE HC + ClinGen HI=3 → 直接 8 分。SVI 指出過度樂觀，應依情境降級：

| 情境 | 判定 | 分數 |
|------|------|------|
| nonsense/frameshift/splice±1,2 且會被 NMD 降解 | `PVS1` | 8 |
| 逃過 NMD + 落在功能域 | `PVS1_Strong` | 4 |
| 逃過 NMD + 移除 >10% 蛋白 | `PVS1_Strong` | 4 |
| 逃過 NMD + 只截尾 | `PVS1_Moderate` | 2 |
| `start_lost`（SVI 上限） | `PVS1_Moderate` | 2 |

- 支撐資料：VEP 掛官方 **NMD plugin**（`NMD.pm` 單檔 bind，不必重建容器）＋ `--total_length`
  （讓 `Protein_position` 變 `123/456` 以算截斷比例）。
- 新增 `PVS1_STRENGTH` / `PVS1_REASON` 欄，可直接篩「被降級的 LoF」。
- **已知簡化**（都寫在 `PVS1_REASON`）：關鍵功能區以 VEP `DOMAINS` 非空代理；未實作「下游 LoF
  在族群中常見」分支；in-frame exon skipping 以 NMD 預測＋截斷比例近似。
- ⚠️ **會改變判讀**：逃過 NMD 又只截尾者可能由 `Likely_Pathogenic` 降為 `VUS`。
- `PVS1_STRENGTH` 有值但 `ACMG_CLASS=Benign` 不是矛盾 —— BA1 是 stand-alone benign 會蓋過一切。

### DGX-2 profile + GPU lock + Pangolin CPU/GPU

- 新增 `dgx` profile：ref_dir/sif_dir 沿用二級 DGX 那份，code 放 `tertiary_code/`，
  掛載 `/datalake_Intermediate,/datalake_Raw,/raid`，`process_high` 48 cores。
- `PANGOLIN_SCORE` 接上二級同一組 `gpu_lock.sh`／`gpu_unlock.sh`（`trap EXIT` 保證還卡），
  `maxForks = 6`（六張 V100）。卡號不寫死，由 lock 動態分配。
- `use_gpu_pangolin` 原本**宣告了卻沒有任何 .nf 使用**（純裝飾）。現已接上：關閉時 export
  空的 `CUDA_VISIBLE_DEVICES` 讓 PyTorch 走 CPU，且 profile 不加 `--nv`，可部署到無 GPU 環境。

### 踩雷：Groovy 布林與 params 覆蓋順序

- **`as boolean` 對字串一律 true**：`--academic_dbnsfp false` 會被判成開啟。一律用
  `params.X.toString().toLowerCase() == 'true'`。三處（snv_annotation / parse_csq / banner）
  必須用同一種寫法，否則 `DBNSFP_VERSION` 會與實際使用的檔案不一致。
- **全域 params 會蓋掉 profile**（本檔全域區塊在 profiles 之後，晚出現者勝）：`use_gpu_lock`
  一度被宣告在全域 `false`，會讓 dgx 的 `true` 失效。**只在 profile 內宣告**。
  （`use_gpu_pangolin` 同樣被覆蓋，但因為沒人使用所以一直沒被發現。）
- **banner 印錯資料庫**：原本寫死印 `params.dbnsfp`，開了 academic flag 仍顯示 4.9c。
  評鑑需要能證明用了哪個資料庫 → 改印實際使用的檔案＋模式。

### 驗證（NA12878 WES，`--academic_dbnsfp true`）

81 欄、`DBNSFP_VERSION=5.3a`；REVEL/MutPred2/VEST4/CADD 各約 25–26%（全表），在 **missense SNV
母體內** 為 REVEL 95.9%、CADD 99.5%、P-KNN 99.5%、MutPred2 98.9%、AlphaMissense 98.8%；
非 missense（synonymous/intron/UTR/indel）一律 `.`，符合 dbNSFP 只收 nsSNV 的設計。
未命中的 missense 多為 PRAMEF 家族等旁系同源區，屬各工具自身覆蓋限制。
ClinGen 對照 77 筆（AGREE 71 / DIFFER_TIER 4 / DIFFER 2）；PVS1 分級 40 `PVS1` + 5 `PVS1_Strong`。
