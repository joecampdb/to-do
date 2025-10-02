# todo

I’m shifting from a **structure-based computational biology** background (proteins, docking, MD) into the **sequence-based -omics** world.  
This document is my **self-study roadmap** — transparent, fast-paced, and pipeline-driven — so people can follow my progress in bioinformatics, especially **genomics & transcriptomics**.

---

## Phase 1: Core NGS (Month 1–2)
- [ ] Learn FASTQ handling & QC  
  - Tools: **FastQC, MultiQC**
- [ ] Practice read trimming & adapter removal  
  - Tools: **Trim Galore, Cutadapt**
- [ ] Master alignment to reference genomes  
  - Tools: **BWA-MEM, Samtools**
- [ ] Post-processing essentials  
  - Tools: **Samtools (sort, index), Picard (mark dups)**
- [ ] Variant calling workflow  
  - Tools: **GATK (HaplotypeCaller, BQSR)**

---

## Phase 2: RNA-Seq (Month 3–4)
- [ ] Bulk RNA quantification  
  - Tools: **Salmon, Kallisto, STAR**
- [ ] Build expression matrices  
- [ ] Differential expression analysis  
  - Tools: **DESeq2, edgeR**
- [ ] Enrichment & pathway analysis  
  - Tools: **clusterProfiler, GSEA**

---

## Phase 3: Single-Cell (Month 5–6)
- [ ] Preprocessing of 10x datasets  
  - Tools: **Cell Ranger**
- [ ] Clustering & visualization  
  - Tools: **Seurat (R), Scanpy (Python)**
- [ ] Learn QC, trajectory inference, and marker detection

---

## Phase 4: Spatial Transcriptomics (Month 7–8)
- [ ] Tissue-mapped RNA-seq  
  - Tools: **Space Ranger (10x), Seurat Spatial**
- [ ] Explore image + transcriptome integration  
  - Tools: **Squidpy**
- [ ] Benchmark datasets from 10x Genomics

---

## Guiding Principles
- 🧬 **Sequence-first**: build intuition around DNA/RNA data.  
- ⚡ **Fast-paced**: 2 months per phase, with checkpoints.  
- 🔄 **Transparent**: share progress openly on social media.  
- 🔗 **Bridge**: connect my structural biology foundation with modern genomics pipelines.  

---

✅ Goal: By ~9 months, be fluent in **NGS → RNA-seq → single-cell → spatial** pipelines and able to contribute to real-world projects in genomics and transcriptomics.
