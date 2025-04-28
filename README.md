# RNA-Seq Analysis Pipeline

This is a **Snakemake**-based RNA-Seq pipeline for performing quality control, read alignment, counting, and differential expression analysis.

---

## 📂 Folder Structure
```
├── adapters/          # Contains adapter sequences for trimming
├── config/            # Contains config.yaml and Samples.tsv
├── raw_data/          # Contains downloaded FASTQ files
├── reference/         # Contains genome reference files (FASTA, GTF, GFF)
├── rules/             # Contains all modular Snakemake rules (*.smk)
├── scripts/           # Contains custom scripts (e.g., DESeq2.R)
└── Snakefile          # Main workflow file
```

---

## 📥 Download Required Files

1. **Reference GFF3 File:**
   - Download manually from:  
     [Homo_sapiens.GRCh38.113.gff3.gz](https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.gff3.gz)
   - Place it inside the `reference/` directory.

2. **Raw FASTQ Files:**
   - Download the raw sequencing data from ENA:  
     [PRJNA229998 - EBI ENA](https://www.ebi.ac.uk/ena/browser/view/PRJNA229998)
   - Place all downloaded `.fastq.gz` files inside the `raw_data/` directory.

---

## ⚙️ Configuration
Edit the `config/config.yaml` file to adjust:
- Quality trimming parameters
- Reference genome paths
- Resource allocations (threads, memory)

Sample metadata (samples and conditions) should be listed in `config/Samples.tsv`.

---

## 🚀 Running the Pipeline

First, perform a dry run:

```bash
snakemake -n -p
```

To run the full pipeline on available cores (example with 8 cores):

```bash
snakemake --cores 8
```

(Optional) Create a full HTML report after completion:

```bash
snakemake --report report.html
```

---

## 🔧 Requirements

- [Snakemake](https://snakemake.readthedocs.io/)
- [Conda](https://docs.conda.io/en/latest/) (recommended)
- Python 3.8+
- R with DESeq2 installed

(You can also set up an `environment.yaml` for reproducibility.)

---

## 📜 Rules Overview

- **quality_control.smk**: Quality and adapter trimming
- **fastqc.smk**: Quality checking using FastQC
- **genome_index.smk**: STAR genome index building
- **star_alignment.smk**: Read alignment to reference genome
- **featurecounts.smk**: Gene expression quantification
- **deseq2_init.smk**: Preparation for DESeq2 analysis
- **deseq2_analysis.smk**: Differential expression analysis
- **visualization.smk**: PCA, MA-plots, Volcano plots

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).

---

## ✨ Acknowledgments

- ENA for raw RNA-Seq data
- Ensembl for the Homo sapiens genome annotations
- Authors of Snakemake, STAR, FastQC, FeatureCounts, and DESeq2
