
rule quality_control:
    input:
        r1 = lambda wildcards: sample_table.loc[wildcards.sample, "fq1"],
        r2 = lambda wildcards: sample_table.loc[wildcards.sample, "fq2"]
    output:
        trimmed1="{sample}/trimmed/{sample}_1.trim.fastq.gz",
        trimmed2="{sample}/trimmed/{sample}_2.trim.fastq.gz",
        unpaired1="{sample}/trimmed/{sample}_1un.trim.fastq.gz",
        unpaired2="{sample}/trimmed/{sample}_2un.trim.fastq.gz",
        summary="{sample}/trimmed/{sample}.txt"
    params:
        qual=config["qual"]["leading"],
        minlen=config["min_len"],
        seed=config["seed_mismatches"],
        palindromeClip=["palindrome_clip_threshold"],
        simpleClip=config["simple_clip_threshold"],
        adapLen=config["adapter_len"],
        adap=config["adapter_file"]
    resources:
        mem_mb=config["resources"]["qc"]["mem_mb"],
        disk_mb=config["resources"]["qc"]["disk_mb"]
    threads: config["resources"]["qc"]["threads"]
    log: "logs/trim/{sample}.log"
    conda: "envs/qc.yaml"
    shell:
        "trimmomatic PE -threads {threads} -summary {output.summary} {input.r1} {input.r2} "
        "{output.trimmed1} {output.unpaired1} {output.trimmed2} {output.unpaired2} "
        "LEADING:{params.qual} TRAILING:{params.qual} MINLEN:{params.minlen} "
        "ILLUMINACLIP:{params.adap}:{params.seed}:{params.palindromeClip}:"
        "{params.simpleClip}:{params.adapLen}:true"
