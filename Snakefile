import pathlib as pl


HG38ACC = "GCF_000001405.40"
HG38FNA = "GCF_000001405.40_GRCh38.p14_genomic.fna"
SAMPLES = [
    f.name.replace(".fastq.gz", "") for f in pl.Path("data/orig/dna_reads").iterdir()
]
RNAs = [
    f.name.replace(".fastq.gz", "") for f in pl.Path("data/orig/rna_reads").iterdir()
]
SAMPLEIDS = set([s[:-2] for s in SAMPLES])


rule all:
    input:
        "data/hg38/README.md",
        "data/hg38/md5sum.txt",
        "data/hg38/ncbi_dataset/fetch.txt",
        "data/hg38/ncbi_dataset/data/dataset_catalog.json",
        "data/hg38/ncbi_dataset/data/assembly_data_report.jsonl",
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/{HG38FNA}",
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/genomic.gff",
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/rna.fna",
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/rna.fna.fai",
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/{HG38FNA}.fai",
        *[f"data/hisat/index/{HG38FNA}/index.{i}.ht2" for i in range(1, 9)],
        *[f"data/hisat/index/rna.fna/index.{i}.ht2" for i in range(1, 9)],
        expand(
            "data/qc/fastqc/orig/dna_reads/{sample}_fastqc.{ext}",
            sample=SAMPLES,
            ext=["zip", "html"],
        ),
        expand(
            "data/qc/fastqc/trim/dna_reads/{sample}/{sample}.{dir}{pair}_fastqc.{ext}",
            sample=SAMPLEIDS,
            dir=["f", "r"],
            pair=["p", "u"],
            ext=["zip", "html"],
        ),
        expand(
            "data/qc/fastqc/trim/rna_reads/{sample}/{sample}.se_fastqc.{ext}",
            sample=RNAs,
            ext=["zip", "html"],
        ),
        expand(
            "data/qc/fastqc/orig/rna_reads/{sample}_fastqc.{ext}",
            sample=RNAs,
            ext=["zip", "html"],
        ),
        expand(
            "data/trim/dna_reads/{sample}/{sample}.{dir}{pair}.fastq.gz",
            sample=SAMPLEIDS,
            dir=["f", "r"],
            pair=["p", "u"],
        ),
        expand(
            "data/trim/rna_reads/{sample}/{sample}.se.fastq.gz",
            sample=RNAs,
        ),
        expand(
            f"data/hisat/mapping/dna_reads/{HG38FNA}/{{sample}}/metrics.txt",
            sample=SAMPLEIDS,
        ),
        expand(
            f"data/hisat/mapping/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.txt",
            sample=SAMPLEIDS,
        ),
        expand(
            f"data/hisat/mapping/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.bam",
            sample=SAMPLEIDS,
        ),
        expand(
            f"data/hisat/mapping/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.bam.bai",
            sample=SAMPLEIDS,
        ),
        expand(
            f"data/hisat/mapping/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.flagstat",
            sample=SAMPLEIDS,
        ),
        expand(
            f"data/hisat/mapping/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.samstats",
            sample=SAMPLEIDS,
        ),
        expand(
            f"data/hisat/mapping/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.bamstats",
            sample=SAMPLEIDS,
        ),
        expand(
            f"data/hisat/mapping/rna_reads/rna.fna/{{sample}}/metrics.txt", sample=RNAs
        ),
        expand(
            f"data/hisat/mapping/rna_reads/rna.fna/{{sample}}/{{sample}}.txt",
            sample=RNAs,
        ),
        expand(
            f"data/hisat/mapping/rna_reads/rna.fna/{{sample}}/{{sample}}.bam",
            sample=RNAs,
        ),
        expand(
            f"data/hisat/mapping/rna_reads/rna.fna/{{sample}}/{{sample}}.bam.bai",
            sample=RNAs,
        ),
        expand(
            f"data/hisat/mapping/rna_reads/rna.fna/{{sample}}/{{sample}}.flagstat",
            sample=RNAs,
        ),
        expand(
            f"data/hisat/mapping/rna_reads/rna.fna/{{sample}}/{{sample}}.samstats",
            sample=RNAs,
        ),
        expand(
            f"data/hisat/mapping/rna_reads/rna.fna/{{sample}}/{{sample}}.bamstats",
            sample=RNAs,
        ),
        expand(
            f"data/hg38.bcf.d/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.bcf",
            sample=SAMPLEIDS,
        ),
        expand(
            f"data/hg38.bcf.d/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.bcfstats",
            sample=SAMPLEIDS,
        ),
        expand(
            f"data/hg38.bcf.d/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.bcf.csi",
            sample=SAMPLEIDS,
        ),
        expand("data/hg38.bcf.d/rna_reads/rna.fna/{sample}/{sample}.bcf", sample=RNAs),
        expand(
            "data/hg38.bcf.d/rna_reads/rna.fna/{sample}/{sample}.bcfstats", sample=RNAs
        ),
        expand(
            "data/hg38.bcf.d/rna_reads/rna.fna/{sample}/{sample}.bcf.csi", sample=RNAs
        ),
        "data/qc/multiqc/multiqc_report.html",


rule calling:
    input:
        ref=f"data/hg38/ncbi_dataset/data/{HG38ACC}/{{fna}}",
        bam="data/hisat/mapping/{mol}/{fna}/{sample}/{sample}.bam",
    output:
        "data/hg38.bcf.d/{mol}/{fna}/{sample}/{sample}.bcf",
    log:
        mpileup="logs/hg38.bcf/{mol}/{fna}/{sample}/mpileup.log",
        callerr="logs/hg38.bcf/{mol}/{fna}/{sample}/calls.log",
        filter="logs/hg38.bcf/{mol}/{fna}/{sample}/filter.log",
    threads: 5
    shell:
        "mkdir -p data/hg38.bcf.d && "
        "bcftools mpileup "
        "--threads {threads} "
        "2> {log.mpileup} "
        "--max-depth 1500 "
        "--fasta-ref {input.ref} "
        "{input.bam} | "
        "bcftools call "
        "--threads {threads} "
        "2> {log.callerr} "
        "--multiallelic-caller "
        "--variants-only | "
        "bcftools filter "
        "2> {log.filter} "
        "--threads {threads} "
        "--include 'QUAL>30 && DP>50' "
        "--output {output} "


rule calls_stats:
    input:
        "data/hg38.bcf.d/{mol}/{fna}/{sample}/{sample}.bcf",
    output:
        "data/hg38.bcf.d/{mol}/{fna}/{sample}/{sample}.bcfstats",
    log:
        "logs/hg38.bcf/{mol}/{fna}/{sample}/stats.log",
    threads: 5
    # need to rename ID for multiqc
    shell:
        "bcftools stats "
        "--threads {threads} "
        "2> {log} "
        "{input} | "
        "sed \"s|{input}|$(echo {input} | cut -d '/' -f 6 | cut -d '.' -f 1)|g\" "
        "> {output} "


rule calls_index:
    input:
        "data/hg38.bcf.d/{mol}/{fna}/{sample}/{sample}.bcf",
    output:
        "data/hg38.bcf.d/{mol}/{fna}/{sample}/{sample}.bcf.csi",
    log:
        "logs/hg38.bcf/{mol}/{fna}/{sample}/index.log",
    threads: 5
    shell:
        "bcftools index "
        "--csi "
        "--threads {threads} "
        "2> {log} "
        "--output {output} "
        "{input}"


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
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/{HG38FNA}",
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/genomic.gff",
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/rna.fna",
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
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/{{ref}}.fna.fai",
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
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/{{fna}}",
    output:
        *[f"data/hisat/index/{{fna}}/index.{i}.ht2" for i in range(1, 9)],
    threads: 8
    resources:
        mem_mb=1024,
    log:
        stdout="logs/hisat/index/{fna}/{fna}.stdout.log",
        stderr="logs/hisat/index/{fna}/{fna}.stderr.log",
    shell:
        "hisat2-build "
        "> {log.stdout} "
        "2> {log.stderr} "
        "--threads {threads} "
        "{input} "
        "data/hisat/index/{wildcards.fna}/index"


rule trimmomatic_pe:
    input:
        fi="data/orig/{mol}/{sample}_1.fastq.gz",
        ri="data/orig/{mol}/{sample}_2.fastq.gz",
    output:
        fpo="data/trim/{mol}/{sample}/{sample}.fp.fastq.gz",
        fuo="data/trim/{mol}/{sample}/{sample}.fu.fastq.gz",
        rpo="data/trim/{mol}/{sample}/{sample}.rp.fastq.gz",
        ruo="data/trim/{mol}/{sample}/{sample}.ru.fastq.gz",
    log:
        stdout="logs/trim/{mol}/{sample}/{sample}.stdout.log",
        stderr="logs/trim/{mol}/{sample}/{sample}.stderr.log",
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


rule trimmomatic_se:
    input:
        "data/orig/{mol}/{sample}.fastq.gz",
    output:
        "data/trim/{mol}/{sample}/{sample}.se.fastq.gz",
    log:
        stdout="logs/trim/{mol}/{sample}/{sample}.stdout.log",
        stderr="logs/trim/{mol}/{sample}/{sample}.stderr.log",
    threads: 5
    resources:
        mem_mb=1024,
    shell:
        "trimmomatic SE "
        "> {log.stdout} "
        "2> {log.stderr} "
        "-phred33 "
        "-threads {threads} "
        "{input} "
        "{output} "
        "TRAILING:20 MINLEN:50"


rule hisat2_refmap:
    input:
        *[f"data/hisat/index/{{fna}}/index.{i}.ht2" for i in range(1, 9)],
        fw="data/trim/{mol}/{sample}/{sample}.fp.fastq.gz",
        rv="data/trim/{mol}/{sample}/{sample}.rp.fastq.gz",
    output:
        met="data/hisat/mapping/{mol}/{fna}/{sample}/metrics.txt",
        sum="data/hisat/mapping/{mol}/{fna}/{sample}/{sample}.txt",
        bam="data/hisat/mapping/{mol}/{fna}/{sample}/{sample}.bam",
        samsorttmp=temp(directory("data/samtools.sort.tmp.d/{mol}.{fna}.{sample}.d/")),
    threads: 20
    resources:
        mem_mb=1024,
    log:
        stderr="logs/hisat/mapping/{mol}/{fna}/{sample}/stderr.log",
        samsort="logs/hisat/mapping/{mol}/{fna}/{sample}/samsort.stderr.log",
    shell:
        "mkdir -p {output.samsorttmp} && "
        "hisat2 "
        "2> {log.stderr} "
        "--fr "
        "--threads {threads} "
        "--time "
        "--no-spliced-alignment "
        "--met-file {output.met} "
        "--summary-file {output.sum} "
        "-x data/hisat/index/{wildcards.fna}/index "
        "-1 {input.fw} "
        "-2 {input.rv} | "
        "samtools sort "
        "-@ {threads} "
        "-O BAM "
        "-o {output.bam} "
        "-T {output.samsorttmp}/tmp "
        "2> {log.samsort}"


rule hisat2_refmap_unpaired:
    input:
        *[f"data/hisat/index/{{fna}}/index.{i}.ht2" for i in range(1, 9)],
        read="data/trim/{mol}/{sample}/{sample}.se.fastq.gz",
    output:
        met="data/hisat/mapping/{mol}/{fna}/{sample}/metrics.txt",
        sum="data/hisat/mapping/{mol}/{fna}/{sample}/{sample}.txt",
        bam="data/hisat/mapping/{mol}/{fna}/{sample}/{sample}.bam",
        samsorttmp=temp(directory("data/samtools.sort.tmp.d/{mol}.{fna}.{sample}.d/")),
    threads: 20
    resources:
        mem_mb=1024,
    log:
        stderr="logs/hisat/mapping/{mol}/{fna}/{sample}/stderr.log",
        samsort="logs/hisat/mapping/{mol}/{fna}/{sample}/samsort.stderr.log",
    shell:
        "mkdir -p {output.samsorttmp} && "
        "hisat2 "
        "2> {log.stderr} "
        "--threads {threads} "
        "--time "
        "--no-spliced-alignment "
        "--met-file {output.met} "
        "--summary-file {output.sum} "
        "-x data/hisat/index/{wildcards.fna}/index "
        "-U {input.read} | "
        "samtools sort "
        "-@ {threads} "
        "-O BAM "
        "-o {output.bam} "
        "-T {output.samsorttmp}/tmp "
        "2> {log.samsort}"


rule hisat2_refmap_samindex:
    input:
        "data/hisat/mapping/{mol}/{fna}/{sample}/{sample}.bam",
    output:
        "data/hisat/mapping/{mol}/{fna}/{sample}/{sample}.bam.bai",
    threads: 5
    log:
        "logs/hisat/mapping/{mol}/{fna}/{sample}/samindex.stderr.log",
    shell:
        "samtools index "
        "2> {log} "
        "-@ {threads} "
        "--output {output} "
        "{input}"


rule hisat2_refmap_samstats:
    input:
        "data/hisat/mapping/{mol}/{fna}/{sample}/{sample}.bam",
    output:
        "data/hisat/mapping/{mol}/{fna}/{sample}/{sample}.samstats",
    threads: 5
    log:
        "logs/hisat/mapping/{mol}/{fna}/{sample}/samstats.stderr.log",
    shell:
        "samtools stats "
        "2> {log} "
        "-@ {threads} "
        "> {output} "
        "{input}"


rule hisat2_refmap_flagstat:
    input:
        "data/hisat/mapping/{mol}/{fna}/{sample}/{sample}.bam",
    output:
        "data/hisat/mapping/{mol}/{fna}/{sample}/{sample}.flagstat",
    threads: 5
    log:
        "logs/hisat/mapping/{mol}/{fna}/{sample}/flagstat.stderr.log",
    shell:
        "samtools flagstat "
        "2> {log} "
        "-@ {threads} "
        "> {output} "
        "{input}"


rule hisat2_refmap_bamstats:
    input:
        "data/hisat/mapping/{mol}/{fna}/{sample}/{sample}.bam",
    output:
        "data/hisat/mapping/{mol}/{fna}/{sample}/{sample}.bamstats",
    threads: 1
    log:
        "logs/hisat/mapping/{mol}/{fna}/{sample}/bamstats.stderr.log",
    shell:
        "bamtools stats "
        "2> {log} "
        "> {output} "
        "-in {input}"


rule fastqc:
    input:
        "data/{src}/{mol}/{sample}.fastq.gz",
    output:
        html="data/qc/fastqc/{src}/{mol}/{sample}_fastqc.html",
        # the suffix _fastqc.zip is necessary for multiqc to find the file
        zip="data/qc/fastqc/{src}/{mol}/{sample}_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "logs/fastqc/{src}/{mol}/{sample}.log",
    threads: 5
    resources:
        mem_mb=1024,
    wrapper:
        "v5.2.1/bio/fastqc"


rule multiqc:
    input:
        expand(
            "data/qc/fastqc/orig/dna_reads/{sample}_fastqc.zip",
            sample=SAMPLES,
        ),
        expand(
            "data/qc/fastqc/trim/dna_reads/{sample}/{sample}.{dir}{pair}_fastqc.zip",
            sample=SAMPLEIDS,
            dir=["f", "r"],
            pair=["p", "u"],
        ),
        expand(
            "data/qc/fastqc/trim/rna_reads/{sample}/{sample}.se_fastqc.zip",
            sample=RNAs,
        ),
        expand(
            "data/qc/fastqc/orig/rna_reads/{sample}_fastqc.zip",
            sample=RNAs,
        ),
        expand(
            "logs/trim/dna_reads/{sample}/{sample}.{std}.log",
            std=["stdout", "stderr"],
            sample=SAMPLEIDS,
        ),
        expand(
            "logs/trim/rna_reads/{sample}/{sample}.{std}.log",
            std=["stdout", "stderr"],
            sample=RNAs,
        ),
        expand(
            f"data/hisat/mapping/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.txt",
            sample=SAMPLEIDS,
        ),
        expand(
            f"data/hisat/mapping/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.flagstat",
            sample=SAMPLEIDS,
        ),
        expand(
            f"data/hisat/mapping/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.samstats",
            sample=SAMPLEIDS,
        ),
        expand(
            f"data/hisat/mapping/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.bamstats",
            sample=SAMPLEIDS,
        ),
        expand(
            f"data/hisat/mapping/rna_reads/rna.fna/{{sample}}/{{sample}}.txt",
            sample=RNAs,
        ),
        expand(
            f"data/hisat/mapping/rna_reads/rna.fna/{{sample}}/{{sample}}.flagstat",
            sample=RNAs,
        ),
        expand(
            f"data/hisat/mapping/rna_reads/rna.fna/{{sample}}/{{sample}}.samstats",
            sample=RNAs,
        ),
        expand(
            f"data/hisat/mapping/rna_reads/rna.fna/{{sample}}/{{sample}}.bamstats",
            sample=RNAs,
        ),
        expand(
            f"data/hg38.bcf.d/dna_reads/{HG38FNA}/{{sample}}/{{sample}}.bcfstats",
            sample=SAMPLEIDS,
        ),
        expand(
            "data/hg38.bcf.d/rna_reads/rna.fna/{sample}/{sample}.bcfstats", sample=RNAs
        ),
    output:
        "data/qc/multiqc/multiqc_report.html",
        directory("data/qc/multiqc/multiqc_data"),
        directory("data/qc/multiqc/multiqc_plots"),
    log:
        stdout="logs/multiqc/stdout.log",
        stderr="logs/multiqc/stderr.log",
    params:
        fqc_dna=expand(
            "data/qc/fastqc/trim/dna_reads/{sample}",
            sample=SAMPLEIDS,
        ),
        fqc_rna=expand(
            "data/qc/fastqc/trim/rna_reads/{sample}",
            sample=RNAs,
        ),
        refmap_dna=expand(
            f"data/hisat/mapping/dna_reads/{HG38FNA}/{{sample}}",
            sample=sorted(SAMPLEIDS),
        ),
        refmap_rna=expand(
            f"data/hisat/mapping/rna_reads/rna.fna/{{sample}}", sample=sorted(RNAs)
        ),
        bcf_dna=expand(
            f"data/hg38.bcf.d/dna_reads/{HG38FNA}/{{sample}}", sample=sorted(SAMPLEIDS)
        ),
        bcf_rna=expand(
            "data/hg38.bcf.d/rna_reads/rna.fna/{sample}", sample=sorted(RNAs)
        ),
    shell:
        "multiqc "
        "> {log.stdout} "
        "2> {log.stderr} "
        "--force "
        "--no-ansi "
        "--export "
        "--no-megaqc-upload "
        "--outdir data/qc/multiqc/ "
        "data/qc/fastqc/orig/dna_reads "
        "data/qc/fastqc/orig/rna_reads "
        "{params.fqc_dna} "
        "{params.fqc_rna} "
        "logs/trim/dna_reads "
        "logs/trim/rna_reads "
        "{params.bcf_dna} "
        "{params.bcf_rna} "
        "{params.refmap_dna} "
        "{params.refmap_rna} "
