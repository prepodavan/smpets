include: "common.smk"


configfile: "snakemake.yml"


HG38ACC = config["reference"]["acc"]
HG38FNA = config["reference"]["filename"]
SAMPLES = list(sorted(config["samples"]["dna"]["names"]))
SAMPLEIDS = list(sorted(config["samples"]["dna"]["ids"]))
RNAs = list(sorted(config["samples"]["rna"]["names"]))
PAIRED = SAMPLEIDS


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
        # dna
        expand(
            "data/orig/dna_reads/{sample}/{sample}_{mate}_fastqc.{ext}",
            sample=SAMPLEIDS,
            mate=["1", "2"],
            ext=["zip", "html"],
        ),
        expand(
            "data/orig/dna_reads/{sample}/{sample}.{dir}{pair}{ext}",
            sample=SAMPLEIDS,
            dir=["f", "r"],
            pair=["p", "u"],
            ext=[
                ".fastq.gz",
                "_fastqc.zip",
                "_fastqc.html",
            ],
        ),
        expand(
            f"data/orig/dna_reads/{{sample}}/{HG38FNA}/{{sample}}.{{tool}}",
            sample=SAMPLEIDS,
            tool=[
                "hisat_metrics.txt",
                "txt",
            ],
        ),
        expand(
            f"data/orig/dna_reads/{{sample}}/{HG38FNA}/{{sample}}{{filter}}.{{tool}}",
            sample=SAMPLEIDS,
            filter=[
                "",
                ".filtered",
            ],
            tool=[
                "bam",
                "bam.bai",
                "flagstat",
                "samstats",
                "bamstats",
                "bcf",
                "bcfstats",
                "vcf",
                "vep.html",
            ],
        ),
        # rna
        expand(
            "data/orig/rna_reads/{sample}/{sample}_fastqc.{ext}",
            sample=RNAs,
            ext=["zip", "html"],
        ),
        expand(
            "data/orig/rna_reads/{sample}/{sample}.se{ext}",
            sample=RNAs,
            ext=[
                ".fastq.gz",
                "_fastqc.zip",
                "_fastqc.html",
            ],
        ),
        expand(
            "data/orig/rna_reads/{sample}/rna.fna/{sample}.{tool}",
            sample=RNAs,
            tool=[
                "hisat_metrics.txt",
                "txt",
            ],
        ),
        expand(
            "data/orig/rna_reads/{sample}/rna.fna/{sample}{filter}.{tool}",
            sample=RNAs,
            filter=[
                "",
                ".filtered",
            ],
            tool=[
                "bam",
                "bam.bai",
                "flagstat",
                "samstats",
                "bamstats",
                "bcf",
                "bcfstats",
                "vcf",
                "vep.html",
            ],
        ),
        "data/qc/multiqc/multiqc_report.html",


rule tabix:
    input:
        f"data/hg38/ncbi_dataset/data/{HG38ACC}/genomic.gff",
    output:
        "data/hg38.genomic.gff.bgz",
        "data/hg38.genomic.gff.bgz.tbi",
    log:
        "logs/tabix.log",
    shell:
        'grep -v "#" {input} | '
        "sort -k1,1 -k4,4n -k5,5n -t$'\\t' | "
        "bgzip -c > {output} 2> {log} && "
        "tabix -p gff {output} 2> {log}"


rule vep:
    input:
        gff="data/hg38.genomic.gff.bgz",
        tbi="data/hg38.genomic.gff.bgz.tbi",
        bcf="{sample_dir}/{fna}/{sample}.bcf",
        fna=f"data/hg38/ncbi_dataset/data/{HG38ACC}/{{fna}}",
    output:
        calls="{sample_dir}/{fna}/{sample}.vcf",
        stats="{sample_dir}/{fna}/{sample}.vep.html",
        vcfin=temp("data/hg38.bcf.d/{sample_dir}/{fna}/{sample}.vcf"),
    log:
        stdout="logs/vep/{sample_dir}/{fna}/{sample}.vcf.stdout.log",
        stderr="logs/vep/{sample_dir}/{fna}/{sample}.vcf.stderr.log",
        convert_input_stdout="logs/vep/{sample_dir}/{fna}/{sample}.vcf.convert.input.stdout.log",
        convert_input_stderr="logs/vep/{sample_dir}/{fna}/{sample}.vcf.convert.input.stderr.log",
    threads: 10
    shell:
        "bcftools convert "
        "> {log.convert_input_stdout} "
        "2> {log.convert_input_stderr} "
        "--threads {threads} "
        "--output-type v "
        "--output {output.vcfin} "
        "{input.bcf} && "
        "vep "
        "> {log.stdout} "
        "2> {log.stderr} "
        "--fork {threads} "
        "--force_overwrite "
        "--everything "
        "--stats_html "
        "--fasta {input.fna} "
        "--gff {input.gff} "
        "--stats_file {output.stats} "
        "--input_file {output.vcfin} "
        "--output_file {output.calls} "


rule calling:
    input:
        ref=f"data/hg38/ncbi_dataset/data/{HG38ACC}/{{fna}}",
        bam="{sample_dir}/{fna}/{sample}.bam",
    output:
        "{sample_dir}/{fna}/{sample}.bcf",
    log:
        mpileup="logs/calling/{sample_dir}/{fna}/{sample}/mpileup.log",
        callerr="logs/calling/{sample_dir}/{fna}/{sample}/calls.log",
        filter="logs/calling/{sample_dir}/{fna}/{sample}/filter.log",
    threads: 5
    params:
        include=calling_filters,
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
        "--include '{params.include}' "
        "--output {output} "


rule calls_stats:
    input:
        "{sample}.bcf",
    output:
        "{sample}.bcfstats",
    log:
        "logs/calling/{sample}/stats.log",
    threads: 5
    # need to rename ID for multiqc
    shell:
        "bcftools stats "
        "--threads {threads} "
        "2> {log} "
        "{input} | "
        "sed 's|{input}|{wildcards.sample}|g' "
        "> {output} "


rule calls_index:
    input:
        "{sample}.bcf",
    output:
        "{sample}.bcf.csi",
    log:
        "logs/calling/{sample}/index.log",
    threads: 5
    shell:
        "bcftools index "
        "--csi "
        "--threads {threads} "
        "2> {log} "
        "--output {output} "
        "{input}"


rule filter_mapping:
    input:
        "{sample}.bam",
    output:
        "{sample}.filtered.bam",
    log:
        "logs/mapping.filtered/{sample}.log",
    params:
        required=branch(
            evaluate(f"{{sample}} in {list(PAIRED)}"),
            then="--require-flags PROPER_PAIR",
            otherwise="",
        ),
    shell:
        "samtools view "
        "2> {log} "
        "--bam "
        "--exclude-flags UNMAP "
        "{params.required} "
        "{input} "
        " > {output}"


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
        fi="data/{sample}_1.fastq.gz",
        ri="data/{sample}_2.fastq.gz",
    output:
        fpo="data/{sample}.fp.fastq.gz",
        fuo="data/{sample}.fu.fastq.gz",
        rpo="data/{sample}.rp.fastq.gz",
        ruo="data/{sample}.ru.fastq.gz",
    log:
        stdout="logs/{sample}.trimmomatic.pe.stdout.log",
        stderr="logs/{sample}.trimmomatic.pe.stderr.log",
    threads: 10
    resources:
        mem_mb=1024,
    params:
        trimmer=trimmomatic_pe_filters,
    shell:
        "trimmomatic PE "
        "> {log.stdout} "
        "2> {log.stderr} "
        "-phred33 "
        "-threads {threads} "
        "{input.fi} {input.ri} "
        "{output.fpo} {output.fuo} "
        "{output.rpo} {output.ruo} "
        "{params.trimmer}"


rule trimmomatic_se:
    input:
        "data/{sample}.fastq.gz",
    output:
        "data/{sample}.se.fastq.gz",
    log:
        stdout="logs/{sample}.trimmomatic.se.stdout.log",
        stderr="logs/{sample}.trimmomatic.se.stderr.log",
    threads: 5
    resources:
        mem_mb=1024,
    params:
        trimmer=trimmomatic_se_filters,
    shell:
        "trimmomatic SE "
        "> {log.stdout} "
        "2> {log.stderr} "
        "-phred33 "
        "-threads {threads} "
        "{input} "
        "{output} "
        "{params.trimmer}"


rule hisat2_refmap:
    input:
        *[f"data/hisat/index/{{fna}}/index.{i}.ht2" for i in range(1, 9)],
        fw="{sample_dir}/{sample}.fp.fastq.gz",
        rv="{sample_dir}/{sample}.rp.fastq.gz",
    output:
        met="{sample_dir}/{fna}/{sample}.hisat_metrics.txt",
        sum="{sample_dir}/{fna}/{sample}.txt",
        bam="{sample_dir}/{fna}/{sample}.bam",
        samsorttmp=temp(
            directory("data/samtools.sort.tmp.d/{sample_dir}/{fna}.{sample}.d/")
        ),
    threads: 20
    resources:
        mem_mb=1024,
    log:
        stderr="logs/hisat/mapping/{sample_dir}/{fna}/{sample}/stderr.log",
        samsort="logs/hisat/mapping/{sample_dir}/{fna}/{sample}/samsort.stderr.log",
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
        read="{sample_dir}/{sample}.se.fastq.gz",
    output:
        met="{sample_dir}/{fna}/{sample}.hisat_metrics.txt",
        sum="{sample_dir}/{fna}/{sample}.txt",
        bam="{sample_dir}/{fna}/{sample}.bam",
        samsorttmp=temp(
            directory("data/samtools.sort.tmp.d/{sample_dir}/{fna}.{sample}.d/")
        ),
    threads: 20
    resources:
        mem_mb=1024,
    log:
        stderr="logs/hisat/mapping/{sample_dir}/{fna}/{sample}/stderr.log",
        samsort="logs/hisat/mapping/{sample_dir}/{fna}/{sample}/samsort.stderr.log",
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
        "data/{sample}.bam",
    output:
        "data/{sample}.bam.bai",
    threads: 5
    log:
        "logs/{sample}/samindex.stderr.log",
    shell:
        "samtools index "
        "2> {log} "
        "-@ {threads} "
        "--output {output} "
        "{input}"


rule hisat2_refmap_samstats:
    input:
        "data/{sample}.bam",
    output:
        "data/{sample}.samstats",
    threads: 5
    log:
        "logs/{sample}/samstats.stderr.log",
    shell:
        "samtools stats "
        "2> {log} "
        "-@ {threads} "
        "> {output} "
        "{input}"


rule hisat2_refmap_flagstat:
    input:
        "data/{sample}.bam",
    output:
        "data/{sample}.flagstat",
    threads: 5
    log:
        "logs/{sample}/flagstat.stderr.log",
    shell:
        "samtools flagstat "
        "2> {log} "
        "-@ {threads} "
        "> {output} "
        "{input}"


rule hisat2_refmap_bamstats:
    input:
        "data/{sample}.bam",
    output:
        "data/{sample}.bamstats",
    threads: 1
    log:
        "logs/{sample}/bamstats.stderr.log",
    shell:
        "bamtools stats "
        "2> {log} "
        "> {output} "
        "-in {input}"


rule split_paired_reads_into_dirs:
    input:
        "{sample_dir}/{sample}_{sfx}.fastq.gz",
    output:
        "{sample_dir}/{sample}/{sample}_{sfx}.fastq.gz",
    log:
        stdout="logs/split/{sample_dir}/{sample}_{sfx}/stdout.log",
        stderr="logs/split/{sample_dir}/{sample}_{sfx}/stderr.log",
    shell:
        "mkdir "
        "> {log.stdout} 2> {log.stderr} "
        "-p $(dirname {output}) && "
        "cp "
        "> {log.stdout} 2> {log.stderr} "
        "{input} {output}"


rule split_unpaired_reads_into_dirs:
    input:
        "{sample_dir}/{sample}.fastq.gz",
    output:
        "{sample_dir}/{sample}/{sample}.fastq.gz",
    log:
        stdout="logs/split/{sample_dir}/{sample}/stdout.log",
        stderr="logs/split/{sample_dir}/{sample}/stderr.log",
    shell:
        "mkdir "
        "> {log.stdout} 2> {log.stderr} "
        "-p $(dirname {output}) && "
        "cp "
        "> {log.stdout} 2> {log.stderr} "
        "{input} {output}"


rule fastqc:
    input:
        "data/{sample}.fastq.gz",
    output:
        html="data/{sample}_fastqc.html",
        # the suffix _fastqc.zip is necessary for multiqc to find the file
        zip="data/{sample}_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "logs/{sample}.log",
    threads: 20
    resources:
        mem_mb=1024,
    wrapper:
        "v5.2.1/bio/fastqc"


rule multiqc:
    input:
        # dna
        expand(
            "data/orig/dna_reads/{sample}/{sample}_{mate}_fastqc.{ext}",
            sample=SAMPLEIDS,
            mate=["1", "2"],
            ext=["zip", "html"],
        ),
        expand(
            "data/orig/dna_reads/{sample}/{sample}.{dir}{pair}{ext}",
            sample=SAMPLEIDS,
            dir=["f", "r"],
            pair=["p", "u"],
            ext=[
                ".fastq.gz",
                "_fastqc.zip",
                "_fastqc.html",
            ],
        ),
        expand(
            f"data/orig/dna_reads/{{sample}}/{HG38FNA}/{{sample}}.{{tool}}",
            sample=SAMPLEIDS,
            tool=[
                "hisat_metrics.txt",
                "txt",
            ],
        ),
        expand(
            f"data/orig/dna_reads/{{sample}}/{HG38FNA}/{{sample}}{{filter}}.{{tool}}",
            sample=SAMPLEIDS,
            filter=[
                "",
                ".filtered",
            ],
            tool=[
                "bam",
                "bam.bai",
                "flagstat",
                "samstats",
                "bamstats",
                "bcf",
                "bcfstats",
                "vcf",
                "vep.html",
            ],
        ),
        # rna
        expand(
            "data/orig/rna_reads/{sample}/{sample}_fastqc.{ext}",
            sample=RNAs,
            ext=["zip", "html"],
        ),
        expand(
            "data/orig/rna_reads/{sample}/{sample}.se{ext}",
            sample=RNAs,
            ext=[
                ".fastq.gz",
                "_fastqc.zip",
                "_fastqc.html",
            ],
        ),
        expand(
            "data/orig/rna_reads/{sample}/rna.fna/{sample}.{tool}",
            sample=RNAs,
            tool=[
                "hisat_metrics.txt",
                "txt",
            ],
        ),
        expand(
            "data/orig/rna_reads/{sample}/rna.fna/{sample}{filter}.{tool}",
            sample=RNAs,
            filter=[
                "",
                ".filtered",
            ],
            tool=[
                "bam",
                "bam.bai",
                "flagstat",
                "samstats",
                "bamstats",
                "bcf",
                "bcfstats",
                "vcf",
                "vep.html",
            ],
        ),
    output:
        "data/qc/multiqc/multiqc_report.html",
        directory("data/qc/multiqc/multiqc_data"),
        directory("data/qc/multiqc/multiqc_plots"),
    log:
        stdout="logs/multiqc/stdout.log",
        stderr="logs/multiqc/stderr.log",
    params:
        dnas=expand("data/orig/dna_reads/{sample}", sample=SAMPLEIDS),
        dna_stats=expand(f"data/orig/dna_reads/{{sample}}/{HG38FNA}", sample=SAMPLEIDS),
        dna_trims=expand("logs/orig/dna_reads/{sample}", sample=SAMPLEIDS),
        rnas=expand("data/orig/rna_reads/{sample}", sample=RNAs),
        rna_stats=expand("data/orig/rna_reads/{sample}/rna.fna", sample=RNAs),
        rna_trims=expand("logs/orig/rna_reads/{sample}", sample=RNAs),
    shell:
        "multiqc "
        "> {log.stdout} "
        "2> {log.stderr} "
        "--force "
        "--no-ansi "
        "--export "
        "--no-megaqc-upload "
        "--outdir data/qc/multiqc/ "
        "{params.dna_trims} "
        "{params.dnas} "
        "{params.dna_stats} "
        "{params.rna_trims} "
        "{params.rnas} "
        "{params.rna_stats} "
