# 三級分析 Pipeline 開發筆記（v3.6）

**負責人：** 林伯昱（p88124019@gs.ncku.edu.tw）
**最後更新：** 2026-08-03

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
├── pangolin_1.0.0.sif
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
    /data/pylin1991/nf-containers/pangolin_1.0.0.sif \
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
# pangolin_1.0.0.sif
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
    /data/pylin1991/nf-containers/{vep_115,pangolin_1.0.0,tertiary_python_1.0.0,annotsv_3.5.10,pharmcat_3.2.0,stellarpgx_graphtyper2.5.1,optitype_1.3.5,samtools_1.23.1}.sif \
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

# 2) Pangolin 能否吃 V100（compute 7.0）
apptainer exec --nv /datalake_Intermediate/pipeline/nextflow_containers/pangolin_1.0.0.sif \
    python3 -c "import torch; print('cuda', torch.cuda.is_available(), torch.cuda.get_device_name(0))"
#   印不出 cuda True → 改用 --use_gpu_pangolin false 走 CPU（只影響速度）

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
