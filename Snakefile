import pathlib as pl


HG38ACC = "GCF_000001405.40"
HG38FNA = "GCF_000001405.40_GRCh38.p14_genomic.fna"
SAMPLES = [
    f.name.replace(".fastq.gz", "") for f in pl.Path("data/orig/dna_reads").iterdir()
]
SAMPLEIDS = set([s[:-2] for s in SAMPLES])
TRIMMED = list()
HISAT_INDEXES = [f"data/hisat/index.{i}.ht2" for i in range(1, 9)]

for sample in SAMPLEIDS:
    TRIMMED.append(f"{sample}.1.forward.paired")
    TRIMMED.append(f"{sample}.1.forward.unpaired")
    TRIMMED.append(f"{sample}.2.reverse.paired")
    TRIMMED.append(f"{sample}.2.reverse.unpaired")


rule all:
    input:
        *[f"data/trim/dna_reads/{tr}.fastq.gz" for tr in TRIMMED],
        "data/hg38/README.md",
        "data/hg38/md5sum.txt",
        "data/hg38/ncbi_dataset/fetch.txt",
        "data/hg38/ncbi_dataset/data/dataset_catalog.json",
        "data/hg38/ncbi_dataset/data/assembly_data_report.jsonl",
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/{HG38FNA}",
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/genomic.gff",
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/rna.fna",
        f"data/faidx/{HG38FNA[:-4]}.fai",
        "data/faidx/rna.fai",
        *HISAT_INDEXES,
        expand(
            "data/qc/fastqc/orig/{sample}_fastqc.{ext}",
            sample=SAMPLES,
            ext=["zip", "html"],
        ),
        expand(
            "data/qc/fastqc/trim/{sample}_fastqc.{ext}",
            sample=TRIMMED,
            ext=["zip", "html"],
        ),
        "data/qc/multiqc/multiqc_report.html",


rule hg38_dehydrated:
    output:
        "data/hg38/README.md",
        "data/hg38/md5sum.txt",
        "data/hg38/ncbi_dataset/fetch.txt",
        "data/hg38/ncbi_dataset/data/dataset_catalog.json",
        "data/hg38/ncbi_dataset/data/assembly_data_report.jsonl",
    log:
        stdout="logs/datasets/hg38dehydrated.stdout.log",
        stderr="logs/datasets/hg38dehydrated.stderr.log",
    threads: 1
    shell:
        "datasets download genome accession "
        f"{HG38ACC} "
        "> {log.stdout} "
        "2> {log.stderr} "
        "--include genome,gff3,rna,seq-report "
        "--reference "
        "--no-progressbar "
        "--filename data/hg38.zip "
        "--dehydrated "
        "&& unzip "
        "> {log.stdout} "
        "2> {log.stderr} "
        "-o -d data/hg38 data/hg38.zip  "
        "&& rm data/hg38.zip"


rule hg38:
    input:
        "data/hg38/ncbi_dataset/fetch.txt",
    output:
        "data/hg38/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna",
        "data/hg38/ncbi_dataset/data/GCF_000001405.40/genomic.gff",
        "data/hg38/ncbi_dataset/data/GCF_000001405.40/rna.fna",
    log:
        stdout="logs/datasets/hg38.stdout.log",
        stderr="logs/datasets/hg38.stderr.log",
    threads: 5
    shell:
        "datasets rehydrate "
        "> {log.stdout} "
        "2> {log.stderr} "
        "--max-workers {threads} "
        "--no-progressbar "
        "--directory data/hg38"


rule faidx:
    input:
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/{{ref}}.fna",
    output:
        "data/faidx/{ref}.fai",
    log:
        stdout="logs/faidx/{ref}.stdout.log",
        stderr="logs/faidx/{ref}.stderr.log",
    threads: 1
    shell:
        "samtools faidx "
        "> {log.stdout} "
        "2> {log.stderr} "
        "--threads {threads} "
        "--fai-idx {output} "
        "{input}"


rule hisat2_build:
    input:
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/{HG38FNA}",
    output:
        *HISAT_INDEXES,
    threads: 8
    resources:
        mem_mb=1024,
    log:
        stdout="logs/hisat/index.stdout.log",
        stderr="logs/hisat/index.stderr.log",
    shell:
        "hisat2-build "
        "> {log.stdout} "
        "2> {log.stderr} "
        "--threads {threads} "
        "{input} "
        "data/hisat/index"


rule trimmomatic_pe:
    input:
        fi="data/orig/dna_reads/{sample}_1.fastq.gz",
        ri="data/orig/dna_reads/{sample}_2.fastq.gz",
    output:
        fpo="data/trim/dna_reads/{sample}.1.forward.paired.fastq.gz",
        fuo="data/trim/dna_reads/{sample}.1.forward.unpaired.fastq.gz",
        rpo="data/trim/dna_reads/{sample}.2.reverse.paired.fastq.gz",
        ruo="data/trim/dna_reads/{sample}.2.reverse.unpaired.fastq.gz",
    log:
        stdout="logs/trim/{sample}.stdout.log",
        stderr="logs/trim/{sample}.stderr.log",
    threads: 10
    resources:
        mem_mb=1024,
    shell:
        "trimmomatic PE "
        "> {log.stdout} "
        "2> {log.stderr} "
        "-phred33 "
        "-threads {threads} "
        "{input.fi} {input.ri} "
        "{output.fpo} {output.fuo} "
        "{output.rpo} {output.ruo} "
        "TRAILING:20 MINLEN:50"


rule fastqc:
    input:
        "data/{src}/dna_reads/{sample}.fastq.gz",
    output:
        html="data/qc/fastqc/{src}/{sample}_fastqc.html",
        # the suffix _fastqc.zip is necessary for multiqc to find the file
        zip="data/qc/fastqc/{src}/{sample}_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "logs/fastqc/{sample}.{src}.log",
    threads: 5
    resources:
        mem_mb=1024,
    wrapper:
        "v5.2.1/bio/fastqc"


rule multiqc:
    input:
        expand("data/qc/fastqc/orig/{sample}_fastqc.zip", sample=SAMPLES),
        expand("data/qc/fastqc/trim/{sample}_fastqc.zip", sample=TRIMMED),
        expand(
            "logs/trim/{sample}.{std}.log", std=["stdout", "stderr"], sample=SAMPLEIDS
        ),
    output:
        "data/qc/multiqc/multiqc_report.html",
        directory("data/qc/multiqc/multiqc_data"),
        directory("data/qc/multiqc/multiqc_plots"),
    log:
        stdout="logs/multiqc/stdout.log",
        stderr="logs/multiqc/stderr.log",
    shell:
        "multiqc "
        "> {log.stdout} "
        "2> {log.stderr} "
        "--force "
        "--no-ansi "
        "--export "
        "--no-megaqc-upload "
        "--outdir data/qc/multiqc/ "
        "data/qc/fastqc/orig/ "
        "data/qc/fastqc/trim/ "
        "logs/trim/ "
