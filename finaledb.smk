import os
import math
import logging

# Snakemake 8 and onwards: "Remote providers have been replaced by Snakemake storage plugins. Please use the corresponding storage plugin instead (snakemake-storage-plugin-*)."
# from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
# HTTP = HTTPRemoteProvider()


def eprint(*args, **kwargs):
    import sys

    sys.stderr.write('\x1b[33m')
    print(*args, file=sys.stderr, **kwargs)
    sys.stderr.write('\x1b[0m')
 
# >>> Configuration >>>
MEM_PER_CORE = int(config.get("MEM_PER_CORE", 1800))

TRIMMOMATIC_CORES = int(config.get("TRIMMOMATIC_CORES", 8))
BWA_CORES = int(config.get("BWA_CORES", 16))
BAM_SORT_CORES = int(config.get("BAM_SORT_CORES", 8))
CRAM_CORES = int(config.get("CRAM_CORES", 4))
INFER_FRAG_CORES = int(config.get("INFER_FRAG_CORES", 8))
INFER_FRAG_SORTBED_CORES = max(1, int(config.get("INFER_FRAG_SORTBED_CORES", INFER_FRAG_CORES - 3)))

WALL_TIME_MAX = int(config.get("WALL_TIME_MAX", 2880))
WALL_TIME_MIN = int(config.get("WALL_TIME_MIN", 300))

WITH_CIGAR = bool(config.get("WITH_CIGAR", 0))
WITH_READ_ID =  bool(config.get("WITH_READ_ID", 0))

# Extra instructions for trimmomatic
TRIMMOMATIC_EXTRA = config.get("TRIMMOMATIC_EXTRA", "")
TRIMMOMATIC_HEAD = config.get("TRIMMOMATIC_HEAD", "")
# <<< Configuration <<<


def get_ref_genome(ref_genome):
    ref_map = {
        "hg19": [os.path.normpath(v) for v in expand(f"data/hg19/human_g1k_v37.fa.gz{{ext}}", ext=["", ".amb", ".ann", ".bwt", ".fai", ".pac", ".sa"])],
        "hg38": [os.path.normpath(v) for v in expand(f"data/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz{{ext}}", ext=["", ".amb", ".ann", ".bwt", ".fai", ".pac", ".sa"])],
    }
    return ref_map[ref_genome]


def get_ref_genome_sizes(ref_genome):
    ref_map = {
        "hg19": os.path.normpath(f"data/hg19/human_g1k_v37.chrom.sizes"),
        "hg38": os.path.normpath(f"data/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set"),
    }
    return ref_map[ref_genome]


# Perform adapter-trimming and filtering using Trimmomatic
rule trimmomatic:
    input: 
        fastq=expand("fastq/{{sample}}.{read_name}.fastq.gz", read_name=["R1", "R2"]),
        adapter="supplementary/TruSeq3-PE-2.fa" # Download from https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE-2.fa"
    output:
        trim=expand("trimmomatic/{{sample}}.R{rg}.{pairing}.fastq.gz", rg=[1, 2], pairing=["paired", "unpaired"]),
    log: "trimmomatic/{sample}.trimmomatic.summary.txt",
    params:
        slurm_job_label=lambda wildcards: f"trimmomatic.{wildcards.sample}",
        trimmer=lambda wildcards, input: config.get("TRIMMER", TRIMMOMATIC_HEAD + f" ILLUMINACLIP:{input.adapter}:2:30:10:2:true TRAILING:5 AVGQUAL:20 MINLEN:36 ") + TRIMMOMATIC_EXTRA,
    threads: lambda wildcards, attempt: int(TRIMMOMATIC_CORES * (0.5 + attempt * 0.5))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        trimmomatic PE -threads {TRIMMOMATIC_CORES} {input.fastq[0]} {input.fastq[1]} $tmpdir/R1.paired.gz $tmpdir/R1.unpaired.gz $tmpdir/R2.paired.gz $tmpdir/R2.unpaired.gz {params.trimmer} 2> >(tee {log} >&2)

        mv $tmpdir/R1.paired.gz {output.trim[0]}.tmp
        mv $tmpdir/R1.unpaired.gz {output.trim[1]}.tmp
        mv $tmpdir/R2.paired.gz {output.trim[2]}.tmp
        mv $tmpdir/R2.unpaired.gz {output.trim[3]}.tmp

        mv {output.trim[0]}.tmp {output.trim[0]}
        mv {output.trim[1]}.tmp {output.trim[1]}
        mv {output.trim[2]}.tmp {output.trim[2]}
        mv {output.trim[3]}.tmp {output.trim[3]}
        """


# Align raw reads using bwa-mem
rule bwa_raw:
    input: 
        ref=lambda wildcards: get_ref_genome(wildcards.ref_genome),
        fastq=expand("trimmomatic/{{sample}}.R{rg}.paired.fastq.gz", rg=[1, 2]),
    output: 
        bam=temp("temp/{sample}.{ref_genome,[^\\.]+}.raw.bam"),
    log:
        log="bam/{sample}.{ref_genome}.samblaster.log",
    params:
        slurm_job_label=lambda wildcards: f"bwa_raw.{wildcards.sample}.{wildcards.ref_genome}",
        extra=lambda wildcards: f"-R '@RG\\tID:{config.get('BWA_RG_ID', wildcards.sample)}\\tSM:{config.get('BWA_SM_ID', wildcards.sample)}\\tPL:{config.get('BWA_PL', 'ILLUMINA')}'"
    threads: lambda wildcards, attempt: int(BWA_CORES * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u
        
        bwa mem -t {BWA_CORES} {params.extra} {input.ref[0]} {input.fastq} | \\
        samblaster 2> >(tee {log} >&2) | samtools view -o $tmpdir/output.bam -
        mv $tmpdir/output.bam {output}
        """


# Sort the BAM file with duplicates marked
rule bam_sort:
    input: 
        bam="temp/{sample}.{ref_genome}.raw.bam",
    output:
        bam="bam/{sample}.{ref_genome,[^\\.]+}.mdups.bam",
        bai="bam/{sample}.{ref_genome}.mdups.bam.bai",
    params:
        slurm_job_label=lambda wildcards: f"bam_sort.{wildcards.sample}.{wildcards.ref_genome}",
        sort_mem=lambda wildcards: int(MEM_PER_CORE * 0.8),
    threads: lambda wildcards, attempt: int(BAM_SORT_CORES * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        samtools sort -@ {BAM_SORT_CORES} -m {params.sort_mem}M -T $tmpdir -o $tmpdir/output.bam {input.bam}
        
        samtools index -@ $(({BAM_SORT_CORES}-1)) $tmpdir/output.bam

        mv $tmpdir/output.bam {output.bam}
        mv $tmpdir/output.bam.bai {output.bai}
        """


rule cram_archive:
    input: 
        ref=lambda wildcards: get_ref_genome(wildcards.ref_genome),
        bam="bam/{sample}.{ref_genome}.mdups.bam",
    output: 
        cram="bam/{sample}.{ref_genome,[^\\.]+}.mdups.cram",
        crai="bam/{sample}.{ref_genome}.mdups.cram.crai",
    params:
        slurm_job_label=lambda wildcards: f"cram_archive.{wildcards.sample}.{wildcards.ref_genome}",
    threads: lambda wildcards, attempt: int(CRAM_CORES * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        samtools view -@ $(({CRAM_CORES}-1)) -T {input.ref[0]} -C -o $tmpdir/output.cram {input.bam}
        samtools index $tmpdir/output.cram

        mv $tmpdir/output.cram {output.cram}
        mv $tmpdir/output.cram.crai {output.crai}
        """


# Infer fragments from query-grouped BAM files
rule infer_fragments:
    input:
        bam="temp/{sample}.{ref_genome}.raw.bam",
        sorted_bam="bam/{sample}.{ref_genome}.mdups.bam",
        sorted_bai="bam/{sample}.{ref_genome}.mdups.bam.bai",
    output:
        frag="frag/{sample}.{ref_genome,[^\\.]+}.frag.bed.gz",
        frag_idx="frag/{sample}.{ref_genome}.frag.bed.gz.tbi",
    log: "frag/{sample}.{ref_genome,[^\\.]+}.log",
    threads: lambda wildcards, attempt: int(INFER_FRAG_CORES * (0.5 + 0.5 * attempt))
    params:
        slurm_job_label=lambda wildcards: f"infer_fragments.{wildcards.sample}.{wildcards.ref_genome}",
        sortbed_mem=lambda wildcards, threads, resources: max(1000, int(INFER_FRAG_SORTBED_CORES * MEM_PER_CORE - 2000)),
        annotate_cmd=lambda wildcards, input: f"annotate_frag.py --bam {input.sorted_bam} | " if WITH_CIGAR else "",
        remove_read_id_cmd=lambda wildcards: "" if WITH_READ_ID else "awk -F'\\t' -v OFS='\\t' '{$4=\".\";print}' | ",
        head_print=lambda wildcards: '"#chrom", "start", "end", "name", "mapq", "strand"',
        head_print2=lambda wildcards: ', "cigar1", "cigar2"' if WITH_CIGAR else "",
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MAX,
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        echo '' | awk -F'\\t' -v OFS='\\t' '{{print {params.head_print}{params.head_print2}}}' | bgzip > $tmpdir/frag.bed.gz
        
        samtools view -h -f 3 -F 3852 -G 48 --incl-flags 48 {input.bam} | \\
        bamToBed -bedpe -mate1 -i stdin | \\
        awk -F'\\t' -v OFS='\\t' '{{if ($1!=$4) next; if ($9=="+") {{s=$2;e=$6}} else {{s=$5;e=$3}} if (e>s) print $1,s,e,$7,$8,$9}}' | \\
        sort -k1,1V -k2,2n --parallel {INFER_FRAG_SORTBED_CORES} -S {params.sortbed_mem}M -T $tmpdir | \\
        {params.annotate_cmd} {params.remove_read_id_cmd} bgzip >> $tmpdir/frag.bed.gz

        tabix -p bed $tmpdir/frag.bed.gz

        mv $tmpdir/frag.bed.gz {output.frag}
        mv $tmpdir/frag.bed.gz.tbi {output.frag_idx}
        """
