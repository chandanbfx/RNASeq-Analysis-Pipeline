
rule star_alignment:
    input:
        trimmed1="{sample}/trimmed/{sample}_1.trim.fastq.gz",
        trimmed2="{sample}/trimmed/{sample}_2.trim.fastq.gz",
        index_tar="genome/index.tar.gz",
        gtf=config["reference"]["gtf"]
    output:
        bam="{sample}/output/{sample}.Aligned.sortedByCoord.out.bam",
        tar_file="{sample}/output.tar.gz",
        check=touch("genome/reference/index/{sample}.check.txt")
    params:
        tar="{sample}/output",
        index="genome/",
        index2="genome/reference/index"
    resources:
        mem_mb=config["resources"]["star"]["mem_mb"],
        disk_mb=config["resources"]["star"]["disk_mb"]
    threads: config["resources"]["star"]["threads"]
    log: "logs/star_alignment/{sample}.log"
    conda: "envs/alignment.yaml"
    shell:
        "tar -xzf {input.index_tar} -C {params.index} && "
        "STAR --runThreadN {threads} --genomeDir {params.index2} "
        "--readFilesIn {input.trimmed1} {input.trimmed2} --readFilesCommand zcat "
        "--sjdbGTFfile {input.gtf} --outFileNamePrefix {params.tar}/{wildcards.sample}. "
        "--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts && "
        "tar -czf {output.tar_file} {params.tar}"
