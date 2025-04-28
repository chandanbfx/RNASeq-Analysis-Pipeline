
rule deseq2_init:
    input:
        counts=expand("quantification/Deseq2Input/{sample}.txt", sample=SAMPLES)
    output:
        "deseq2/all.rds"
    params:
        directory="quantification/Deseq2Input/",
        sample_table = config["sample_table"]
    resources:
        mem_mb=config["resources"]["deseq2"]["mem_mb"],
        disk_mb=config["resources"]["deseq2"]["disk_mb"]
    threads: config["resources"]["deseq2"]["threads"]
    conda:
        "envs/deseq2.yaml"    
    log: "logs/deseq2/init.log"
    script: "scripts/deseq2Init.R"
