
rule fastqc:
    input:
        trimmed1="{sample}/trimmed/{sample}_1.trim.fastq.gz",
        trimmed2="{sample}/trimmed/{sample}_2.trim.fastq.gz"
    output:
        fastqc_1="{sample}/trimmed/{sample}_1.trim_fastqc.html",
        fastqc_2="{sample}/trimmed/{sample}_2.trim_fastqc.html",
        fastqc_zip1="{sample}/trimmed/{sample}_1.trim_fastqc.zip",
        fastqc_zip2="{sample}/trimmed/{sample}_2.trim_fastqc.zip"
    resources:
        mem_mb=config["resources"]["qc"]["mem_mb"],
        disk_mb=config["resources"]["qc"]["disk_mb"]
    threads: config["resources"]["qc"]["threads"]
    log: "logs/fastqc/{sample}.log"
    conda:
         "envs/qc.yaml"
    shell:
        "fastqc -t {threads} {input.trimmed1} {input.trimmed2} --dir /data1"
