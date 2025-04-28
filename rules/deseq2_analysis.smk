

rule deseq2:
    input:
        "deseq2/all.rds"
    output:
        table="deseq2/results/diffexp/{contrast}.diffexp.tsv",
        ma_plot="deseq2/results/diffexp/{contrast}.ma-plot.svg"
    params:
        contrast = lambda wildcards: wildcards.contrast
    resources:
        mem_mb=config["resources"]["deseq2"]["mem_mb"],
        disk_mb=config["resources"]["deseq2"]["disk_mb"]
    threads: config["resources"]["deseq2"]["threads"]
    conda:
        "envs/deseq2.yaml"
    log: "logs/deseq2/{contrast}.diffexp.log"
    script: "scripts/DESeq2.R"
