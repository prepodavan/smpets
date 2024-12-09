import pathlib as pl


SAMPLES = [
    f.name.replace(".fastq.gz", "") for f in pl.Path("data/orig/dna_reads").iterdir()
]
SAMPLEIDS = set([s[:-2] for s in SAMPLES])
TRIMMED = list()

for sample in SAMPLEIDS:
    TRIMMED.append(f"{sample}.1.forward.paired")
    TRIMMED.append(f"{sample}.1.forward.unpaired")
    TRIMMED.append(f"{sample}.2.reverse.paired")
    TRIMMED.append(f"{sample}.2.reverse.unpaired")


rule all:
    input:
        *[f"data/trim/dna_reads/{tr}.fastq.gz" for tr in TRIMMED],
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
