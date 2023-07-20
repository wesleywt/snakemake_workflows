configfile: "config.yaml"

def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule all:
    input:
        "results/plots/quals.svg"

rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        temp("results/mapped_reads/{sample}.bam")
    params:
        rg=r'@RG\tID:{sample}\tSM:{sample}'
    log:
        "logs/bwe_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb -> {output}) 2> {log}"


rule samtools_sort:
    input:
        "results/mapped_reads/{sample}.bam"
    output:
        protected("results/sorted_reads/{sample}.bam")
    shell:
        "samtools sort -T results/sorted_reads/{wildcards.sample} -O bam {input} > {output}"

rule samtools_index:
    input:
        "results/sorted_reads/{sample}.bam"
    output:
        "results/sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("results/sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("results/sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "results/calls/all_vcf"
    params:
        prm = config["prior_mutation_rate"]
    log:
        "logs/calls/all_vcf.log"
    shell:
        "(bcftools mpileup -P {params.prm} -f {input.fa} {input.bam} | bcftools call -mv -> {output}) 2> {log}"

rule plot_quals:
    input:
        'results/calls/all_vcf'
    output:
        'results/plots/quals.svg'
    script:
        'scripts/plot-quals.py'
