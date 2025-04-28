
rule genome_index:
    input:
        fa=config["reference"]["genome"],
        gtf=config["reference"]["gtf"]
    output:
        index_tar="genome/index.tar.gz"
    params:
        tar="genome/index",
        chrBin=config["star"]["chr_bin_nbits"],
        sjdbOverhang=config["star"]["sjdb_overhang"]
    resources:
        mem_mb=config["resources"]["star"]["mem_mb"],
        disk_mb=config["resources"]["star"]["disk_mb"]
    threads: config["resources"]["star"]["threads"]
    log: "logs/star_index/index.log"
    conda: "envs/alignment.yaml"
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate "
        "--genomeDir {params.tar} --genomeFastaFiles {input.fa} "
        "--sjdbGTFfile {input.gtf} --genomeChrBinNbits={params.chrBin} "
        "--sjdbOverhang={params.sjdbOverhang} && "
        "tar -czf {output.index_tar} {params.tar}"
