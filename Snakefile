# Snakefile: RNA-seq pipeline from SRA accessions to DE and enrichment results
import csv
import os

with open("metadata/samples.csv") as f:
    SAMPLES = [row["sample"] for row in csv.DictReader(f)]

GENCODE_RELEASE = "50"
TX_URL = (f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
          f"release_{GENCODE_RELEASE}/gencode.v{GENCODE_RELEASE}.transcripts.fa.gz")


rule all:
    input:
        "results/qc/multiqc_report.html",
        "results/sessionInfo.txt"


rule download_fastq:
    output:
        r1 = temp("raw_data/{sample}_1.fastq"),
        r2 = temp("raw_data/{sample}_2.fastq")
    threads: 4
    log: "logs/download/{sample}.log"
    conda: "envs/sra.yaml"
    shell:
        "fasterq-dump {wildcards.sample} --split-files --threads {threads} "
        "--outdir raw_data > {log} 2>&1"


rule fastqc:
    input:
        r1 = "raw_data/{sample}_1.fastq",
        r2 = "raw_data/{sample}_2.fastq"
    output:
        html1 = "results/qc/fastqc/{sample}_1_fastqc.html",
        zip1  = "results/qc/fastqc/{sample}_1_fastqc.zip",
        html2 = "results/qc/fastqc/{sample}_2_fastqc.html",
        zip2  = "results/qc/fastqc/{sample}_2_fastqc.zip"
    threads: 2
    log: "logs/fastqc/{sample}.log"
    conda: "envs/qc.yaml"
    shell:
        "fastqc -t {threads} -o results/qc/fastqc {input.r1} {input.r2} > {log} 2>&1"


rule download_transcriptome:
    output: "ref/gencode_transcripts.fa.gz"
    params: url = TX_URL
    log: "logs/download_transcriptome.log"
    conda: "envs/download.yaml"
    shell:
        "curl -L -o -f {output} {params.url} > {log} 2>&1"


rule salmon_index:
    input: "ref/gencode_transcripts.fa.gz"
    output: directory("ref/salmon_index")
    threads: 8
    log: "logs/salmon_index.log"
    conda: "envs/salmon.yaml"
    shell:
        "salmon index -t {input} -i {output} -p {threads} > {log} 2>&1"


rule salmon_quant:
    input:
        r1 = "raw_data/{sample}_1.fastq",
        r2 = "raw_data/{sample}_2.fastq",
        index = "ref/salmon_index"
    output:
        "quants/{sample}_quant/quant.sf"
    params:
        outdir=lambda w, output: os.path.dirname(output[0])
    threads: 4
    log: "logs/salmon/{sample}.log"
    conda: "envs/salmon.yaml"
    shell:
        "salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} "
        "-p {threads} -o {params.outdir} > {log} 2>&1"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_{read}_fastqc.zip",
               sample=SAMPLES, read=["1", "2"]),
        expand("quants/{sample}_quant/quant.sf", sample=SAMPLES)
    output:
        "results/qc/multiqc_report.html"
    log: "logs/multiqc.log"
    conda: "envs/qc.yaml"
    shell:
        "multiqc -f -o results/qc -n multiqc_report.html "
        "results/qc/fastqc quants > {log} 2>&1"


rule r_analysis:
    input:
        quants=expand("quants/{sample}_quant/quant.sf",sample=SAMPLES),
        samples="metadata/samples.csv",
        scripts=expand("scripts/{s}",s=["config.R", "04_deseq2.R",
                                        "05_qc_plots.R",
                                        "06_enrichment.R", "run_all.R"])
    output:
        "results/tables/DE_amlodipine_vs_control.csv",
        "results/sessionInfo.txt"
    log: "logs/r_analysis.log"
    conda: "envs/r.yaml"
    shell:
        "Rscript scripts/run_all.R > {log} 2>&1"