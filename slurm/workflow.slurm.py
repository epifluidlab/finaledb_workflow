# FinaleDB pipeline

__author__ = "Haizi Zheng"
__copyright__ = "Copyright 2021, Haizi Zheng"
__email__ = "haizi.zh@gmail.com"
__license__ = "MIT"

from os import path
import os
import sys


# smk file home location
smk_path = path.normpath(path.dirname(__file__))
sys.path = [smk_path] + sys.path

from os import path
from tempfile import TemporaryDirectory
from snakemake.shell import shell
from snakemake.logging import logger


def test():
    shell("echo $(which fastqc) > {snakemake.output}")


def fastqc():
    # input: expand("data/xenograft-20201201/fastq/{{fq_id}}_{rg}.fq.gz", rg=[1, 2])
    # output: expand("results/20201130-week-49/xenograft/fastqc/{{fq_id}}_{rg}_fastqc.{ext}", rg=[1, 2], ext=["html", "zip"])
    with TemporaryDirectory() as tempdir:
        job_label = snakemake.params.get("label")
        if job_label:
            logger.info(f"Job started: {job_label}")
        logger.info(f"Using temporary directory: {tempdir}")

        cmd = f"fastqc -t {snakemake.threads} -o {tempdir} -f fastq {snakemake.input}"
        logger.info(cmd)
        shell(cmd)

        destdir = os.path.dirname(snakemake.output[0])
        logger.info(f"Copying to destination: {destdir}")

        for f in os.listdir(tempdir):
            if f.endswith("fastqc.html") or f.endswith("fastqc.zip"):
                logger.info("Copying {f}")
                full_path = path.join(tempdir, f)
                shell("mv {full_path} {destdir}")


fastqc_post_trim = fastqc


def trimmomatic():
    trimmomatic_args = snakemake.config.get(
        "TRIMMOMATIC",
        (
            f"HEADCROP:6 ILLUMINACLIP:{snakemake.input.adapter}:2:30:10:2:keepBothReads "
            "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:36 "
        ),
    )
    logger.info("Trimmomatic arguments:")
    logger.info(trimmomatic_args)

    with TemporaryDirectory() as tempdir:
        job_label = snakemake.params.get("label")
        if job_label:
            logger.info(f"Job started: {job_label}")
        logger.info(f"Using temporary directory: {tempdir}")

        cmd = (
            f"trimmomatic PE -threads {snakemake.threads} "
            f"{snakemake.input.fastq[0]} {snakemake.input.fastq[1]} "
            f"{tempdir}/R1P.fq.gz {tempdir}/R1U.fq.gz "
            f"{tempdir}/R2P.fq.gz {tempdir}/R2U.fq.gz "
            f"{trimmomatic_args} "
            f"2>&1 | tee {snakemake.output.log}"
        )
        logger.info(cmd)
        shell(cmd)

        logger.info(f"Copying to destination")
        shell("mv {tempdir}/R1P.fq.gz {snakemake.output.trim[0]}.tmp")
        shell("mv {tempdir}/R1U.fq.gz {snakemake.output.trim[1]}.tmp")
        shell("mv {tempdir}/R2P.fq.gz {snakemake.output.trim[2]}.tmp")
        shell("mv {tempdir}/R2U.fq.gz {snakemake.output.trim[3]}.tmp")

        shell("mv {snakemake.output.trim[0]}.tmp {snakemake.output.trim[0]}")
        shell("mv {snakemake.output.trim[1]}.tmp {snakemake.output.trim[1]}")
        shell("mv {snakemake.output.trim[2]}.tmp {snakemake.output.trim[2]}")
        shell("mv {snakemake.output.trim[3]}.tmp {snakemake.output.trim[3]}")


def bwa_raw():
    with TemporaryDirectory() as tempdir:
        job_label = snakemake.params.get("label")
        if job_label:
            logger.info(f"Job started: {job_label}")
        logger.info(f"Using temporary directory: {tempdir}")

        sample = snakemake.wildcards.sample
        rgid = sample
        rgpu = rgid
        rglb = sample

        bwa_threads = max(snakemake.threads - 2, 1)
        logger.info(f"Total cores available: {snakemake.threads}")
        logger.info(f"# of cores for bwa-mem: {bwa_threads}")

        cmd = (
            f"bwa mem -t {bwa_threads} "
            f"-R '@RG\\tID:{rgid}\\tSM:{sample}\\tPL:ILLUMINA\\tLB:{rglb}\\tPU:{rgpu}' "
            f"{snakemake.input.ref[0]} {snakemake.input.fastq} | "
            f"samblaster 2> >(tee {tempdir}/samblaster.log >&2) | "
            f"samtools view -b - > "
            f"{tempdir}/output.bam "
        )
        logger.info(cmd)

        shell(cmd)

        logger.info(f"Copying to destination")
        shell("mv {tempdir}/output.bam {snakemake.output.bam}.tmp")
        shell("mv {snakemake.output.bam}.tmp {snakemake.output.bam}")
        shell("mv {tempdir}/samblaster.log {snakemake.output.log}")


def infer_fragments():
    with TemporaryDirectory() as tempdir:
        job_label = snakemake.params.get("label")
        if job_label:
            logger.info(f"Job started: {job_label}")
        logger.info(f"Using temporary directory: {tempdir}")

        nthreads = snakemake.threads
        total_mem = snakemake.resources.mem_mb
        mem_sortbed = int((total_mem - 2000) * 0.8)
        logger.info(f"# of cores available: {nthreads}")
        logger.info(f"Total memory available: {total_mem}MB")
        logger.info(f"Memory used by sort-bed: {mem_sortbed}MB")

        # logger.info(
        #     "Infer fragments using: "
        #     "samtools view -h -f 3 -F 3852 -q 10 | grep -v -e 'XA:Z:' -e 'SA:Z:'"
        # )

        if snakemake.params.get("with_read_id", False):
            perl_script = """perl -ne 'chomp;@f=split " ";if($f[0] ne $f[3]){{next;}}$s=$f[1];$e=$f[5];if($f[8] eq "-"){{$s=$f[4];$e=$f[2];}}if($e>$s){{print "$f[0]\\t$s\\t$e\\t.\\t$f[7]\\t$f[8]\\n";}}' | """
        else:
            perl_script = """perl -ne 'chomp;@f=split " ";if($f[0] ne $f[3]){{next;}}$s=$f[1];$e=$f[5];if($f[8] eq "-"){{$s=$f[4];$e=$f[2];}}if($e>$s){{print "$f[0]\\t$s\\t$e\\t$f[6]\\t$f[7]\\t$f[8]\\n";}}' | """

        shell(
            "samtools view -h -f 3 -F 3852 {snakemake.input.bam} | "
            # "grep -v -e 'XA:Z:' -e 'SA:Z:' | "
            f"bamToBed -bedpe -mate1 -i stdin | {perl_script}"
            "sort-bed --max-mem {mem_sortbed}M --tmpdir {tempdir} - | "
            "bgzip > {tempdir}/output.bed.gz"
        )

        logger.info("Indexing fragments")
        shell("tabix -p bed {tempdir}/output.bed.gz")

        logger.info(f"Copying to destination")
        shell("mv {tempdir}/output.bed.gz {snakemake.output.frag}.tmp")
        shell("mv {snakemake.output.frag}.tmp {snakemake.output.frag}")
        shell("mv {tempdir}/output.bed.gz.tbi {snakemake.output.frag_idx}")


def bam_sort():
    # input:
    #     bam="results/20201130-week-49/xenograft/temp/{sample}.raw.bam",
    # output:
    #     bam="results/20201130-week-49/xenograft/bam/{sample}.mdups.bam",
    #     bai="results/20201130-week-49/xenograft/bam/{sample}.mdups.bam.bai"
    with TemporaryDirectory() as tempdir:
        job_label = snakemake.params.get("label")
        if job_label:
            logger.info(f"Job started: {job_label}")
        logger.info(f"Using temporary directory: {tempdir}")

        nthreads = snakemake.threads
        total_mem = snakemake.resources.mem_mb
        mem_per_thread = int((total_mem - 2000) / nthreads)
        logger.info(f"# of cores available: {nthreads}")
        logger.info(f"Total memory available: {total_mem}MB")
        logger.info(f"Memory per thread allocated to samtools sort: {mem_per_thread}MB")

        shell(
            "samtools sort -@ {snakemake.threads} "
            "-m {mem_per_thread}M -o {tempdir}/output.bam {snakemake.input.bam}"
        )
        logger.info("Indexing BAM file")
        shell("samtools index -@ {snakemake.threads} {tempdir}/output.bam")

        logger.info(f"Copying to destination")
        shell("mv {tempdir}/output.bam {snakemake.output.bam}.tmp")
        shell("mv {snakemake.output.bam}.tmp {snakemake.output.bam}")
        shell("mv {tempdir}/output.bam.bai {snakemake.output.bai}")


def bam_stats():
    with TemporaryDirectory() as tempdir:
        job_label = snakemake.params.get("label")
        if job_label:
            logger.info(f"Job started: {job_label}")
        logger.info(f"Using temporary directory: {tempdir}")

        nthreads = snakemake.threads
        total_mem = snakemake.resources.mem_mb
        logger.info(f"# of cores available: {nthreads}")
        logger.info(f"Total memory available: {total_mem}MB")

        shell("samtools stats {snakemake.input.bam} > {tempdir}/stats.txt")
        shell("samtools idxstats {snakemake.input.bam} > {tempdir}/idxstats.txt")
        shell("samtools flagstats {snakemake.input.bam} > {tempdir}/flagstats.txt")

        logger.info(f"Copying to destination")
        shell("mv {tempdir}/stats.txt {snakemake.output.stats}")
        shell("mv {tempdir}/idxstats.txt {snakemake.output.idxstats}")
        shell("mv {tempdir}/flagstats.txt {snakemake.output.flagstats}")


def picard_metrics_insert_size():
    # Run fastqc, since there can be race conditions if multiple jobs
    # use the same fastqc dir, we create a temp dir.
    with TemporaryDirectory() as tempdir:
        job_label = snakemake.params.get("label")
        if job_label:
            logger.info(f"Job started: {job_label}")
        logger.info(f"Using temporary directory: {tempdir}")

        memory = ""
        if "mem_mb" in snakemake.resources.keys():
            memory = "-Xmx{}M".format(str(snakemake.resources["mem_mb"]))

        temp_txt = path.join(tempdir, "insert_size.txt")
        temp_pdf = path.join(tempdir, "insert_size.pdf")
        temp_log = path.join(tempdir, "insert_size.log")

        java_opts = snakemake.params.get("java_opts", "")
        java_opts = f"{java_opts} -Dpicard.useLegacyParser=false"
        if "-Xmx" not in java_opts:
            java_opts = f"{java_opts} {memory}"

        cmd = (
            f"gatk --java-options '{java_opts}' "
            "CollectInsertSizeMetrics "  # Tool and its subcommand
            f"-I {snakemake.input.bam} "  # Input file
            f"-O {temp_txt} "  # Output metrics
            f"-H {temp_pdf} "  # Output metrics
            f"2> >(tee {temp_log} >&2)"
        )
        logger.info(cmd)
        shell(cmd)

        # Move outputs into proper position.
        logger.info(f"Copying to destination")
        shell("mv {temp_txt} {snakemake.output.insert_size}")
        shell("mv {temp_pdf} {snakemake.output.insert_size_pdf}")
        shell("mv {temp_log} {snakemake.output.insert_size_log}")


def picard_gc_bias():
    # Run fastqc, since there can be race conditions if multiple jobs
    # use the same fastqc dir, we create a temp dir.
    with TemporaryDirectory() as tempdir:
        job_label = snakemake.params.get("label")
        if job_label:
            logger.info(f"Job started: {job_label}")
        logger.info(f"Using temporary directory: {tempdir}")

        memory = ""
        if "mem_mb" in snakemake.resources.keys():
            memory = "-Xmx{}M".format(str(snakemake.resources["mem_mb"]))

        temp_txt = path.join(tempdir, "gc_bias.txt")
        temp_chart = path.join(tempdir, "gc_bias.pdf")
        temp_summary = path.join(tempdir, "gc_bias.summary.txt")
        temp_log = path.join(tempdir, "gc_bias.log")

        java_opts = snakemake.params.get("java_opts", "")
        java_opts = f"{java_opts} -Dpicard.useLegacyParser=false"
        if "-Xmx" not in java_opts:
            java_opts = f"{java_opts} {memory}"

        cmd = (
            f"gatk --java-options '{java_opts}' "
            "CollectGcBiasMetrics "  # Tool and its subcommand
            f"-I {snakemake.input.bam} "  # Input file
            f"-R {snakemake.input.ref} "  # ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
            f"-O {temp_txt} "  # Output metrics
            f"-CHART {temp_chart} "  # Output metrics
            f"-S {temp_summary} "  # Output metrics
            f"2> >(tee {temp_log} >&2)"
        )
        logger.info(cmd)
        shell(cmd)

        # Move outputs into proper position.
        shell("mv {temp_txt} {snakemake.output.gc_bias}")
        shell("mv {temp_chart} {snakemake.output.gc_bias_chart}")
        shell("mv {temp_summary} {snakemake.output.gc_bias_summary}")
        shell("mv {temp_log} {snakemake.output.gc_bias_log}")


globals()[snakemake.rule]()