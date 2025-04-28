
rule quantification:
    input:
        bam="{sample}/output/{sample}.Aligned.sortedByCoord.out.bam",
        gff=config["reference"]["gff"]
    output:
        read_count="quantification/output/{sample}.readCount.txt",
        extracted="quantification/Deseq2Input/{sample}.txt"
    resources:
        mem_mb=config["resources"]["count"]["mem_mb"],
        disk_mb=config["resources"]["count"]["disk_mb"]
    threads: config["resources"]["count"]["threads"]
    log: "logs/featureCount/{sample}.log"
    conda:
         "envs/quantification.yaml"
    shell:
        "featureCounts -T {threads} -t gene -g ID -p -a {input.gff} "
        "-o {output.read_count} {input.bam} && "
        "awk -v OFS='\t' '! /#/{{print $1, $7}}' {output.read_count} | "
        "sed '1d' > {output.extracted}"
