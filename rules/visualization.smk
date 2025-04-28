
rule pca:
    input:
        "deseq2/all.rds"
    output:
        "deseq2/results/pca.svg"
    params:
        pca_labels=config["pca"]["labels"]
    resources:
        mem_mb=config["resources"]["deseq2"]["mem_mb"],
        disk_mb=config["resources"]["deseq2"]["disk_mb"]
    threads: config["resources"]["deseq2"]["threads"]
    conda: "envs/visualization.yaml"
    log: "logs/pca.log"
    script: "scripts/plotPCA.R"

rule heatmap:
    input:
        "deseq2/all.rds"
    output:
        "deseq2/results/heatmap.svg"
    resources:
        mem_mb=config["resources"]["deseq2"]["mem_mb"],
        disk_mb=config["resources"]["deseq2"]["disk_mb"]
    threads: config["resources"]["deseq2"]["threads"]
    conda: "envs/visualization.yaml"
    log: "logs/heatmap.log"
    script: "scripts/plotHeatmap.R"

rule volcanoPlot:
    input:
        "deseq2/results/diffexp/treated-vs-untreated.diffexp.tsv"
    output:
        vol_plot="deseq2/results/diffexp/volcano.svg"
    resources:
        mem_mb=config["resources"]["deseq2"]["mem_mb"],
        disk_mb=config["resources"]["deseq2"]["disk_mb"]
    threads: config["resources"]["deseq2"]["threads"]
    conda: "envs/visualization.yaml"
    log: "logs/deseq2/volcano.log"
    script: "scripts/volcano.R"
