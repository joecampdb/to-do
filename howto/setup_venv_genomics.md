# 🧬 Chapter I — Setting Up an NGS Environment (`omics`)

A complete, step-by-step guide for aspiring bioinformaticians to create a clean environment for Next-Generation Sequencing (NGS) analysis and prepare a human reference genome (GRCh38).

---

## 📦 1. Install Mambaforge

Mambaforge is a lightweight Conda distribution that includes the faster **Mamba** solver.

### For Apple Silicon (M1/M2/M3)
```bash
curl -L -o miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash miniforge.sh -b -p $HOME/mambaforge
source ~/mambaforge/etc/profile.d/conda.sh
```

### For Intel Macs
```bash
curl -L -o miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh
bash miniforge.sh -b -p $HOME/mambaforge
source ~/mambaforge/etc/profile.d/conda.sh
```

### Verify installation
```bash
conda --version
```

If `mamba` is missing:
```bash
conda install -n base -c conda-forge mamba -y
```

---

## ⚙️ 2. Configure Conda Channels

Add channels for bioinformatics tools:
```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

---

## 🧫 3. Create the NGS Environment (`omics`)

Create an environment with essential NGS tools:
```bash
mamba create -n omics   fastqc multiqc cutadapt trim-galore bwa samtools picard gatk4 sra-tools pigz -y
```

Activate it:
```bash
conda activate omics
```

---

## 🧩 4. Install Common Extensions

Inside the `(omics)` environment:
```bash
conda install -c bioconda -c conda-forge   bedtools bcftools vcftools entrez-direct enaBrowserTools seqtk tabix -y
```

Verify key tools:
```bash
fastqc --version
multiqc --version
bwa 2>&1 | head -n 1
samtools --version | head -n 1
picard MarkDuplicates --version
gatk --version
bedtools --version
```

---

## 📁 5. Setup Directory Structure

Organize your workspace:
```bash
mkdir -p ~/omics_project/{data/raw_fastq,ref,results,logs}
cd ~/omics_project
```

Your structure:
```
omics_project/
├── data/
│   └── raw_fastq/
├── ref/
├── results/
└── logs/
```

---

## 🧬 6. Download a Human Reference Genome (GRCh38)

Fetch the reference genome (FASTA) and annotation (GTF):

```bash
cd ~/omics_project/ref

# FASTA (reference genome)
wget -c https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa ref.fa

# GTF (gene annotation)
wget -c https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
gunzip Homo_sapiens.GRCh38.112.gtf.gz
mv Homo_sapiens.GRCh38.112.gtf annotation.gtf
```

---

## 🧠 7. Index the Reference Genome

Build an index for alignment:
```bash
cd ~/omics_project/ref
bwa index ref.fa
samtools faidx ref.fa
```

---

## 🧫 8. (Optional) Test the Environment with a Small Dataset

Download a small E. coli test set:
```bash
cd ~/omics_project/data/raw_fastq

wget -c https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/008/SRR390728/SRR390728.fastq.gz
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
mv GCF_000005845.2_ASM584v2_genomic.fna ref.fa

bwa index ref.fa
bwa mem ref.fa SRR390728.fastq.gz | samtools view -S -b | samtools sort -o aln.bam
samtools index aln.bam
```

Run QC:
```bash
fastqc SRR390728.fastq.gz
multiqc .
```

---

## 🧭 9. Fetching Real Datasets

### From NCBI SRA
```bash
prefetch SRR12345678
fasterq-dump SRR12345678 --split-files -e 8
gzip SRR12345678_*.fastq
```

### From ENA (faster HTTPS)
```bash
enaDataGet -f fastq SRR12345678
```

### Search for Neuroscience Data
```bash
esearch -db sra -query "Homo sapiens[orgn] AND brain[title] AND RNA-Seq" | efetch -format runinfo | head
```

---

## 🧹 10. Deactivate Environment

When finished:
```bash
conda deactivate
```

List all environments:
```bash
conda env list
```

Remove if necessary:
```bash
conda remove -n omics --all
```

---

## ✅ Summary

You now have:
- A clean Conda/Mamba environment (`omics`)
- All core NGS tools (FastQC → GATK)
- Common utilities (bedtools, bcftools, entrez-direct)
- Organized directories for data and results
- A human reference genome (GRCh38) ready for alignment

This completes **Chapter I — Environment Setup for NGS Bioinformatics**.
