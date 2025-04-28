# Snakefile
configfile: "config/config.yaml"

import pandas as pd
from snakemake.utils import validate

# Load samples
sample_table = pd.read_table(config["sample_table"], index_col="sample")
SAMPLES = sample_table.index.tolist()

# Load contrasts
contrasts = list(config["diffexp"]["contrasts"].keys())

# Make available to all rules
def get_sample_table():
    return sample_table

def get_samples():
    return SAMPLES

def get_contrasts():
    return contrasts

print("Loaded samples:", SAMPLES)
print("Loaded contrasts:", contrasts)

conda:
    "envs/base.yaml"

# Include all sub-workflows
include: "rules/quality_control.smk"
include: "rules/fastqc.smk"
include: "rules/genome_index.smk"
include: "rules/star_alignment.smk"
include: "rules/featurecounts.smk"
include: "rules/deseq2_analysis.smk"
include: "rules/deseq2_init.smk"
include: "rules/visualization.smk"

# Final targets
rule all:
    input:
        # Preprocessing outputs
        expand("{sample}/trimmed/{sample}_1.trim_fastqc.html", sample=get_samples()),
        
        # Differential expression outputs
        "deseq2/all.rds",
        expand("deseq2/results/diffexp/{contrast}.diffexp.tsv", contrast=get_contrasts()),
        expand("deseq2/results/diffexp/{contrast}.ma-plot.svg", contrast=get_contrasts()),
        "deseq2/results/pca.svg",
        "deseq2/results/heatmap.svg",
        "deseq2/results/diffexp/volcano.svg"
