# This snakemake specification file defines the workflow working for SLURM environment.

import os
import math

# Some rules, such as bwa, can make the full use of all CPU cores available. By
# default, each node contains 28 cores, which is the RM partition at PSC. Users
# may override this parameter from command line by: --config FULL_CORE=[# of
# cores]
FULL_CORES = config.get("FULL_CORES", 28)

# Available memory (GB) with each CPU core allocated. At PSC, each core comes
# approximately 4.5GB memory.
MEM_MB_PER_CORE = config.get("MEM_MB_PER_CORE", 4200)

# RG ID for bwa mem
BWA_RG_ID = config.get("BWA_ID", "anonymous")

# We perform all computation tasks in Docker containers
docker_image = "docker://zephyre/comp-bio:v0.3.7"

# Path to reference genome assemlies
ref_genome_fa = {
    "hg19": {
        "path": "data/ref_genome/human_g1k_v37.fasta.gz"
    },
    "hg38": {
        "path": "data/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
    }
}


chrom_sizes = {
    "hg38": "data/chrom_sizes/GRCh38_EBV.chrom.sizes",
    "hg19": "data/chrom_sizes/human_g1k_v37.chrom.sizes",
    "ucsc_hg19": "data/chrom_sizes/hg19.chrom.sizes"
}


def nlogn(size_bytes):
    """
    A handy helper function that non-linearly increases as the file size grows larger.
    We use this function to model the maximal time needed while making certain SLURM requests.
    """
    import math

    size_gb = max(size_bytes / 1024**3, 1)
    return (math.log10(size_gb) + 1) * size_gb


def mem_mb_per_core(wildcards, threads):
    """Claim certain amount of memory based on the number of threads"""
    return threads * MEM_MB_PER_CORE


rule trim:
    singularity: docker_image
    input:
        fastq=expand("fastq/{{entry_id}}.{read_name}.fastq.gz", read_name=["R1", "R2"]),
        adapter="data/adapters/TruSeq3-PE-2.fa"
    output:
        fastq=temp(expand("trim/{{entry_id}}.{read_name}.{paired}.fastq.gz", read_name=["R1", "R2"], paired=["paired", "unpaired"])),
        summary="trim/{entry_id}.trimmomatic.summary.txt"
    threads: 3
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=mem_mb_per_core,
        time_min=lambda wildcards, input, threads: round(nlogn(os.path.getsize(input.fastq[0])) * 60 / threads + 30),
        input_size=lambda wildcards, input: round(os.path.getsize(input.fastq[0]) / 1024**2)
    params:
        label=lambda wildcards: f"{wildcards.entry_id}",
        partition="RM-shared"
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local

        echo `date +"%F %T %Z"` "Performing trimmomatic..."
        trimmomatic PE -threads {threads} {input.fastq[0]} {input.fastq[1]} $LOCAL/R1P.fq.gz $LOCAL/R1U.fq.gz $LOCAL/R2P.fq.gz $LOCAL/R2U.fq.gz \
            CROP:50 ILLUMINACLIP:{input.adapter}:2:30:10:2:keepBothReads MINLEN:36 2>&1 | tee $LOCAL/summary.txt

        echo `date +"%F %T %Z"` "Copying rule results to target directory..."
        cp -a $LOCAL/summary.txt {output.summary}

        set +e; rand_str=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 16); set -e
        temp_file="{output.fastq[0]}.$rand_str"
        cp -a "$LOCAL/R1P.fq.gz" "$temp_file"
        mv "$temp_file" {output.fastq[0]}
        set +e; rand_str=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 16); set -e
        temp_file="{output.fastq[1]}.$rand_str"
        cp -a "$LOCAL/R1U.fq.gz" "$temp_file"
        mv "$temp_file" {output.fastq[1]}
        set +e; rand_str=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 16); set -e
        temp_file="{output.fastq[2]}.$rand_str"
        cp -a "$LOCAL/R2P.fq.gz" "$temp_file" 
        mv "$temp_file" {output.fastq[2]}
        set +e; rand_str=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 16); set -e
        temp_file="{output.fastq[3]}.$rand_str"
        cp -a "$LOCAL/R2U.fq.gz" "$temp_file"
        mv "$temp_file" {output.fastq[3]}
        """


rule fastqc:
    singularity: docker_image
    input: expand("trim/{{entry_id}}.R{read_id}.paired.fastq.gz", read_id=[1, 2])
    output: expand("fastqc/{{entry_id,EE[0-9]+}}.R{read_id}.paired_fastqc.{ext}", read_id=[1, 2], ext=["html", "zip"])
    threads: 2
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=mem_mb_per_core,
        time_min=lambda wildcards, input, threads: round(nlogn(os.path.getsize(input[0])) * 20 / threads + 30),
        input_size=lambda wildcards, input: round(os.path.getsize(input[0]) / 1024**2)
    params:
        label=lambda wildcards: f"{wildcards.entry_id}",
        partition="RM-shared"
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local

        echo `date +"%F %T %Z"` "FASTQC started"
        fastqc -t {threads} -o $LOCAL -f fastq {input}
        echo `date +"%F %T %Z"` "FASTQC finished"

        echo `date +"%F %T %Z"` "Copying rule results to target directory..."
        cp -a $LOCAL/{wildcards.entry_id}.R1.paired_fastqc.html {output[0]}
        cp -a $LOCAL/{wildcards.entry_id}.R1.paired_fastqc.zip {output[1]}
        cp -a $LOCAL/{wildcards.entry_id}.R2.paired_fastqc.html {output[2]}
        cp -a $LOCAL/{wildcards.entry_id}.R2.paired_fastqc.zip {output[3]}
        """


rule bwa_raw:
    singularity: docker_image
    input:
        ref=lambda wildcards: ref_genome_fa[wildcards.assembly]["path"],
        ref_index=lambda wildcards: [f"{ref_genome_fa[wildcards.assembly]['path']}.{ext}" \
            for ext in ["fai", "gzi", "sa", "ann", "amb", "pac", "bwt"]],
        fastq=ancient(expand("trim/{{entry_id}}.R{read_id}.paired.fastq.gz", read_id=[1, 2]))
    output:
        bam=temp("temp/{entry_id,EE[0-9]+}.{assembly,[a-zA-Z1-9]+}.raw.bam")
    threads: FULL_CORES
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=mem_mb_per_core,
        time_min=lambda wildcards, input, threads: round(nlogn(os.path.getsize(input.fastq[0])) * 900 / threads + 30),
        input_size=lambda wildcards, input: round(os.path.getsize(input.fastq[0]) / 1024**2)
    params:
        label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}",
        partition="RM"
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local

        echo `date +"%F %T %Z"` "bwa mem mapping..."
        bwa mem -t {threads} -R \'@RG\\tID:{BWA_RG_ID}\\tSM:{wildcards.entry_id}\' {input.ref} {input.fastq} | 
            samtools view -b - > $LOCAL/temp.raw.bam

        set +e; rand_str=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 16); set -e
        temp_file="{output.bam}.$rand_str"
        cp -a "$LOCAL/temp.raw.bam" "$temp_file"
        mv "$temp_file" {output.bam}
        """


# Query-sort the raw BAM file for samblaster dedup
rule dedup_query_sort:
    singularity: docker_image
    input:
        bam="temp/{entry_id,EE[0-9]+}.{assembly,[a-zA-Z1-9]+}.raw.bam"
    output:
        bam=temp("temp/{entry_id,EE[0-9]+}.{assembly,[a-zA-Z1-9]+}.qsorted.bam")
    threads: 4
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4200,
        time_min=lambda wildcards, input, threads: math.ceil(os.path.getsize(input.bam) / 1024**3 * 15 / threads) + 45,
        input_size=lambda wildcards, input: round(os.path.getsize(input.bam) / 1024**2)
    params:
        label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}",
        partition="RM-shared",
        mem_mb_per_thread=lambda wildcards, threads, resources: int(resources.mem_mb * 0.8 / threads)
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local

        samtools sort -@ {threads} -m {params.mem_mb_per_thread}M -n -o $LOCAL/temp.qsorted.bam -T /local/ {input.bam}

        set +e; rand_str=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 16); set -e
        temp_file="{output.bam}.$rand_str"
        cp -a "$LOCAL/temp.qsorted.bam" "$temp_file"
        mv "$temp_file" "{output.bam}"
        """


# Query-sort the raw BAM file for samblaster dedup
rule dedup_samblaster:
    singularity: docker_image
    input:
        bam="temp/{entry_id}.{assembly}.qsorted.bam"
    output:
        bam="bam/{entry_id,EE[0-9]+}.{assembly,[a-zA-Z1-9]+}.mdups.bam",
        bai="bam/{entry_id}.{assembly}.mdups.bam.bai",
        summary="bam/{entry_id}.{assembly}.samblaster.summary.txt"
    threads: 4
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4200,
        time_min=lambda wildcards, input, threads: math.ceil(os.path.getsize(input.bam) / 1024**3 * 30 / threads) + 45,
        input_size=lambda wildcards, input: round(os.path.getsize(input.bam) / 1024**2)
    params:
        label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}",
        partition="RM-shared",
        mem_mb_per_thread=lambda wildcards, threads, resources: int(resources.mem_mb * 0.8 / threads)
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local


        echo `date +"%F %T %Z"` "samblaster deduplicating..."
        samtools view -@ {threads} -h {input.bam} | samblaster 2> >(tee $LOCAL/summary.txt >&2) > $LOCAL/temp.mdups.unsorted.sam

        echo `date +"%F %T %Z"` "sorting..."
        samtools sort -@ {threads} -m {params.mem_mb_per_thread}M -o $LOCAL/output.mdups.bam -T /local/ $LOCAL/temp.mdups.unsorted.sam

        samtools index -@ {threads} $LOCAL/output.mdups.bam $LOCAL/output.mdups.bam.bai

        set +e; rand_str=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 16); set -e
        temp_file="{output.bam}.$rand_str"
        cp -a "$LOCAL/output.mdups.bam" "$temp_file"
        mv "$temp_file" "{output.bam}"

        cp -a "$LOCAL/output.mdups.bam.bai" {output.bai}
        cp -a "$LOCAL/summary.txt" {output.summary}
        """