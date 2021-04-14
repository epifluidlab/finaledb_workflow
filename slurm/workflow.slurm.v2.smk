import os
import math

def eprint(*args, **kwargs):
    import sys

    sys.stderr.write('\x1b[33m')
    print(*args, file=sys.stderr, **kwargs)
    sys.stderr.write('\x1b[0m')
 
 # >>> Configuration >>>
FULL_CORES = config.get("FULL_CORES", 28)
PART_CORES = config.get("PART_CORES", 3)
# <<< Configuration <<<

# # Load the entryId / sraId mapping
# def load_entry_mapping():
#     import csv
#     with open("finaledb.entries.csv", "r") as f:
#         reader = csv.reader(f, delimiter=",")
#         return {f"EE{v[0].strip()}": v[2].strip() for v in reader if v[-1].strip() == "cristiano2019nature"}

# entry_mapping = load_entry_mapping()

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

def format_time_str(total_minutes):
    hour = total_minutes / 60
    minutes = total_minutes % 60
    return "%02d:%02d:00" % (hour, minutes)

def get_file_size(file_path):
    import os

    return os.path.getsize(file_path)


def nlogn(size_bytes):
    import math

    size_gb = max(size_bytes / 1024**3, 1)
    return (math.log10(size_gb) + 1) * size_gb

# def get_unprocessed_bam(entry_id, assembly):
#     sra_id = entry_sra_mapping[entry_id]
#     return f"unprocessed2/{sra_id}/{sra_id}.{assembly}.mdups.bam.gz"

# # Some BAM files are not at the right place
# rule get_bam:
#     resources:
#         disk_mb=100000
#     input:
#         bam=lambda wildcards: get_unprocessed_bam(wildcards.sample_id, wildcards.assembly)
#     output:
#         bam="entries/{sample_id}/{assembly}/{sample_id}.{assembly}.mdups.bam",
#         bai="entries/{sample_id}/{assembly}/{sample_id}.{assembly}.mdups.bam.bai"
#     shell:
#         """
#         gzip -dc {input.bam} > {output.bam}
#         samtools index {output.bam}
#         """


rule bam2fq:
    singularity: "docker://zephyre/comp-bio:v0.3.7"
    input:
        "adalsteinsson/{entry_id}.hg19.mdups.bam"
    output:
        fastq=expand("fastq/{{entry_id,EE[0-9]+}}.{read_name}.fastq.gz", read_name=["R1", "R2"])
    threads: 3
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4200,
        time_min=lambda wildcards, input, threads: round(nlogn(os.path.getsize(input[0])) * 40 / threads + 30),
        input_size=lambda wildcards, input: round(os.path.getsize(input[0]) / 1024**2)
    params:
        mem_mb_per_thread=lambda wildcards, threads, resources: int(resources.mem_mb * 0.8 / threads),
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local

        echo `date +"%F %T %Z"` "Preparing the workspace..."
        cp -a {input} $LOCAL/input.bam

        echo `date +"%F %T %Z"` "[samtools sort] Sorting by query name..."
        samtools sort -@ {threads} -m {params.mem_mb_per_thread}M -n -o $LOCAL/temp.qsort.bam -T $LOCAL $LOCAL/input.bam

        echo `date +"%F %T %Z"` "[bedtools bamtofastq] Extracting reads..."
        mkfifo $LOCAL/temp.R1.fastq.pipe $LOCAL/temp.R2.fastq.pipe
        gzip < $LOCAL/temp.R1.fastq.pipe > $LOCAL/output.R1.fastq.gz &
        gzip < $LOCAL/temp.R2.fastq.pipe > $LOCAL/output.R2.fastq.gz &
        bedtools bamtofastq -i $LOCAL/temp.qsort.bam \
            -fq $LOCAL/temp.R1.fastq.pipe \
            -fq2 $LOCAL/temp.R2.fastq.pipe 
        sleep 20

        echo `date +"%F %T %Z"` "Copying rule results to target directory..."
        cp -a $LOCAL/output.R1.fastq.gz {output.fastq[0]}
        cp -a $LOCAL/output.R2.fastq.gz {output.fastq[1]}
        """

rule trim:
    singularity: "docker://zephyre/comp-bio:v0.3.7"
    input:
        fastq=expand("fastq/{{entry_id}}.{read_name}.fastq.gz", read_name=["R1", "R2"]),
        # fastq=expand("fastq/{{entry_id}}.{read_name}.fastq.bz2", read_name=["R1", "R2"]),
        adapter="data/adapters/TruSeq3-PE-2.fa"
    output:
        fastq=expand("trim/{{entry_id}}.{read_name}.{paired}.fastq.gz", read_name=["R1", "R2"], paired=["paired", "unpaired"]),
        summary="trim/{entry_id}.trimmomatic.summary.txt"
    threads: 3
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4200,
        time_min=lambda wildcards, input, threads: round(nlogn(os.path.getsize(input.fastq[0])) * 60 / threads + 30),
        input_size=lambda wildcards, input: round(os.path.getsize(input.fastq[0]) / 1024**2)
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local

        echo `date +"%F %T %Z"` "Decompressing fastq files..."
        # lbzcat -n {threads} {input.fastq[0]} > $LOCAL/R1.fq
        # lbzcat -n {threads} {input.fastq[1]} > $LOCAL/R2.fq
        pigz -p {threads} -dc {input.fastq[0]} > $LOCAL/R1.fq
        pigz -p {threads} -dc {input.fastq[1]} > $LOCAL/R2.fq

        echo `date +"%F %T %Z"` "Performing trimmomatic..."
        trimmomatic PE -threads {threads} $LOCAL/R1.fq $LOCAL/R2.fq $LOCAL/R1P.fq.gz $LOCAL/R1U.fq.gz $LOCAL/R2P.fq.gz $LOCAL/R2U.fq.gz \
            CROP:50 ILLUMINACLIP:{input.adapter}:2:30:10:2:keepBothReads MINLEN:36 2>&1 | tee $LOCAL/summary.txt

        echo `date +"%F %T %Z"` "Copying rule results to target directory..."
        cp -a $LOCAL/summary.txt {output.summary}
        cp -a $LOCAL/R1P.fq.gz {output.fastq[0]}
        cp -a $LOCAL/R1U.fq.gz {output.fastq[1]}
        cp -a $LOCAL/R2P.fq.gz {output.fastq[2]}
        cp -a $LOCAL/R2U.fq.gz {output.fastq[3]}
        """

rule fastqc:
    singularity: "docker://zephyre/comp-bio:v0.3.7"
    input: expand("trim/{{entry_id}}.R{read_id}.paired.fastq.gz", read_id=[1, 2])
    output: expand("fastqc/{{entry_id,EE[0-9]+}}.R{read_id}.paired_fastqc.{ext}", read_id=[1, 2], ext=["html", "zip"])
    threads: 2
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4200,
        time_min=lambda wildcards, input, threads: round(nlogn(os.path.getsize(input[0])) * 20 / threads + 30),
        input_size=lambda wildcards, input: round(os.path.getsize(input[0]) / 1024**2)
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


# rule bwa_dc_split:
#     # Since bwa mem may take significant time to complete, one way around is splitting the FASTQ files into pieces,
#     # dc means divide-and-conquer
#     # process each of them and combine the 
#     singularity: "docker://zephyre/comp-bio:v0.3.7"
#     input: "trim/{entry_id}.{read_name}.paired.fastq.gz"
#     output:
#         # A manifest of all splitted FASTQ parts. File name example: EEnnnnn.R1.manifest.100m.txt
#         # Chunk size: 100m means each file contains 100m reads
#         manifest="temp/bwa_dc/{entry_id,EE[0-9]+}.{read_name,R[12]}.manifest.{chunk_size}.txt",
#         # split_fastq="temp/bwa_pre_split/{entry_id}.{read_name}.split.{chunk_id,[0-9]+}.fastq"
#     threads: 2
#     resources:
#         cpus=lambda wildcards, threads: threads,
#         mem_mb=lambda wildcards, threads: threads * 4200,
#         time_min=lambda wildcards, input, threads: round(nlogn(os.path.getsize(input[0])) * 8 + 15),
#         input_size=lambda wildcards, input: round(os.path.getsize(input[0]) / 1024**2)
#     params:
#         label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.read_name}.{wildcards.chunk_size}-{wildcards.chunk_id}",
#         partition="RM-shared"
#     shell:
#         """
#         IFS=$'\\n\\t '
#         LOCAL=/local

#         project_home=`pwd`

#         # Determine the chunk size
#         split_size=$((`echo {wildcards.chunk_size} | sed 's/.$//'` * 4000000))
#         echo "Split size: $split_size lines"

#         input_fastq={input[0]}
#         echo `date +"%F %T %Z"` "Splitting $input_fastq ..."
#         input_fastq=`readlink -f $project_home/$input_fastq`
#         ls -lah $input_fastq

#         pushd $LOCAL

#         prefix="{wildcards.entry_id}.{wildcards.read_name}.split.{wildcards.chunk_size}."
#         gzip -dc $input_fastq | split -a 4 --numeric-suffixes=1 -l $split_size - $prefix

#         # Determine the total number of chunks
#         chunk_counter=0

#         part_file_list="$prefix*"
#         for part_file in `echo $part_file_list`; do
#             chunk_counter=$(($chunk_counter + 1))
            
#             part_file=$(basename $part_file)
#             set +e
#             rand_str=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 16)
#             set -e
#             temp_file="$project_home/temp/bwa_dc/${{prefix}}tmp.$rand_str"
#             cp -a $part_file $temp_file
#             mv $temp_file "$project_home/temp/bwa_dc/$part_file.fastq"
#         done

#         echo $chunk_counter > "$project_home/{output.manifest}"
#         """


# rule bwa_dc_map:
#     # Run bwa mem for each splitted chunk
#     singularity: "docker://zephyre/comp-bio:v0.3.7"
#     input: 
#         ref=lambda wildcards: ref_genome_fa[wildcards.assembly]["path"],
#         ref_index=lambda wildcards: [f"{ref_genome_fa[wildcards.assembly]['path']}.{ext}" \
#             for ext in ["fai", "gzi", "sa", "ann", "amb", "pac", "bwt"]],
#         fastq=expand("temp/bwa_dc/{{entry_id}}.{read_name}.split.{{chunk_size}}.{{chunk_id}}.fastq", read_name=["R1", "R2"])
#     output: "temp/bwa_dc/{entry_id}.{assembly}.split.{chunk_size}.{chunk_id,[0-9]+}.sam"
#     threads: FULL_CORES
#     resources:
#         cpus=lambda wildcards, threads: threads,
#         mem_mb=lambda wildcards, threads: threads * 4200,
#         # Time needed: for 1G R1 reads, it requires ~250 SUs
#         time_min=lambda wildcards, input, threads: round(int(wildcards.chunk_size[:-1]) * 1e-3 * 300 * 60 / threads + 30),
#         # time_min=lambda wildcards, input, threads: round(nlogn(os.path.getsize(input.fastq[0])) * 900 / threads + 30),
#         input_size=lambda wildcards, input: round(os.path.getsize(input.fastq[0]) / 1024**2)
#     params:
#         cores_bwa=lambda wildcards, threads, resources: threads,
#         label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}.{wildcards.chunk_size}-{wildcards.chunk_id}",
#         partition="RM"
#     shell:
#         """
#         IFS=$'\\n\\t '
#         LOCAL=/local

#         project_home=`pwd`
#         r1=$(readlink -f {input.fastq[0]})
#         r2=$(readlink -f {input.fastq[1]})

#         pushd $LOCAL

#         bwa mem -t {threads} -R "@RG\\tID:haizi\\tSM:{wildcards.entry_id}" "$project_home/{input.ref}" $r1 $r2 > tmp.sam

#         set +e
#         rand_str=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 16)
#         set -e
#         temp_file="$project_home/{output}.$rand_str"
#         cp -a tmp.sam $temp_file
#         mv $temp_file "$project_home/{output}"
#         """


# # Read the manifest file and return a list of BAM files listed
# def load_manifest(wildcards, chunk_size=None):
#     if not chunk_size:
#         chunk_size = wildcards.chunk_size
#     manifest_files = [f"temp/bwa_dc/{wildcards.entry_id}.{read_name}.manifest.{chunk_size}.txt" for read_name in ["R1", "R2"]]
#     part_nos = []
#     for idx in [0, 1]:
#         with open(manifest_files[idx], "rt") as f:
#             part_nos.append(int(f.read()))

#     # Sanity check: the number of chunks for R1 and R2 should be the same
#     assert part_nos[0] == part_nos[1]

#     return [f"temp/bwa_dc/{wildcards.entry_id}.{wildcards.assembly}.split.{chunk_size}.{str(idx + 1).zfill(4)}.sam" for idx in range(part_nos[0])]


# rule bwa_dc_merge:
#     # Collect all splitted BAM files according to the manifest, and merge to the raw BAM file as a whole
#     singularity: "docker://zephyre/comp-bio:v0.3.7"
#     input:
#         manifest=expand("temp/bwa_dc/{{entry_id}}.{read_name}.manifest.4m.txt", read_name=["R1", "R2"]),
#         sam=lambda wildcards: load_manifest(wildcards, "4m")
#     threads: 4
#     resources:
#         cpus=lambda wildcards, threads: threads,
#         mem_mb=lambda wildcards, threads: threads * 4200,
#         # Time: for 4 cores, each sam file with 1M chunk size needs 1 minute
#         time_min=lambda wildcards, input, threads: int(len(input.sam) * 16 / threads + 10),
#         input_size=lambda wildcards, input: len(input.sam)
#     output:
#         bam=temp("temp/bwa_dc/{entry_id,EE[0-9]+}.{assembly,[a-zA-Z1-9]+}.raw.bam")
#     params:
#         label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}",
#         partition="RM-shared"
#     shell:
#         """
#         IFS=$'\\n\\t '
#         LOCAL=/local

#         set -x
#         samtools merge -@ {threads} -c -p $LOCAL/temp.raw.bam {input.sam}

#         set +e
#         rand_str=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 16)
#         set -e
#         temp_file="{output.bam}.$rand_str"
#         cp -a $LOCAL/temp.raw.bam $temp_file
#         mv $temp_file {output.bam}

#         echo `date +"%F %T %Z"` "Cleaning up..."
#         rm -f {input.sam}
#         """

# # rule bwa_raw:
# #     singularity: "docker://zephyre/comp-bio:v0.3.7"
# #     input:
# #         ref=lambda wildcards: ref_genome_fa[wildcards.assembly]["path"],
# #         ref_index=lambda wildcards: [f"{ref_genome_fa[wildcards.assembly]['path']}.{ext}" \
# #             for ext in ["fai", "gzi", "sa", "ann", "amb", "pac", "bwt"]],
# #         fastq=ancient(expand("trim/{{entry_id}}.R{read_id}.paired.fastq.gz", read_id=[1, 2]))
# #     output:
# #         bam=temp("temp/{entry_id,EE[0-9]+}.{assembly,[a-zA-Z1-9]+}.raw.bam")
# #     threads: FULL_CORES
# #     resources:
# #         cpus=lambda wildcards, threads: threads,
# #         mem_mb=lambda wildcards, threads: threads * 4200,
# #         time_min=lambda wildcards, input, threads: round(nlogn(os.path.getsize(input.fastq[0])) * 900 / threads + 30),
# #         # time_min=lambda wildcards, input, threads: math.ceil(os.path.getsize(input.fastq[0]) / 1024**3 * 900 / threads) + 90,
# #         input_size=lambda wildcards, input: round(os.path.getsize(input.fastq[0]) / 1024**2)
# #     params:
# #         cores_bwa=lambda wildcards, threads, resources: threads
# #     shell:
# #         """
# #         echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
# #         LOCAL=/local

# #         ls -lah {input.fastq[0]}
# #         ls -lah {input.fastq[1]}

# #         echo `date +"%F %T %Z"` "bwa mem mapping..."
# #         bwa mem -t {threads} -R \'@RG\\tID:haizi\\tSM:{wildcards.entry_id}\' {input.ref} {input.fastq} > $LOCAL/temp.raw.sam
# #         #    samtools view -b - > $LOCAL/temp.raw.bam

# #         echo `date +"%F %T %Z"` "Copying rule results to target directory..."
# #         samtools view -@ {threads} -b $LOCAL/temp.raw.sam > {output.bam}

# #         # cp -a $LOCAL/temp.raw.bam {output.bam}
# #         """

rule bwa_samblaster:
    singularity: "docker://zephyre/comp-bio:v0.3.7"
    input:
        bam="temp/bwa_dc/{entry_id,EE[0-9]+}.{assembly,[a-zA-Z1-9]+}.raw.bam"
    output:
        bam="bam/{entry_id,EE[0-9]+}.{assembly,[a-zA-Z1-9]+}.mdups.bam",
        bai="bam/{entry_id}.{assembly,[a-zA-Z1-9]+}.mdups.bam.bai",
        summary="bam/{entry_id}.{assembly}.samblaster.summary.txt"
    threads: 4
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4200,
        time_min=lambda wildcards, input, threads: math.ceil(os.path.getsize(input.bam) / 1024**3 * 30 / threads) + 45,
        input_size=lambda wildcards, input: round(os.path.getsize(input.bam) / 1024**2)
    params:
        mem_mb_per_thread=lambda wildcards, threads, resources: int(resources.mem_mb * 0.8 / threads),
        label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}",
        partition="RM-shared"
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local
        output_basename=$(basename {output.bam})

        echo `date +"%F %T %Z"` "Preparing the workspace"
        ls -lsh {input.bam}
        cp -a {input.bam} $LOCAL/raw.bam

        echo `date +"%F %T %Z"` "queyr-sorting..."

        temp_file="${{output_basename}}.temp.qsorted.bam"
        if [[ -f "temp/bwa_dc/$temp_file" ]]; then
            echo "Found cached file: temp/bwa_dc/$temp_file"
            cp -av "temp/bwa_dc/$temp_file" $LOCAL/temp.qsorted.bam
        else
            samtools sort -@ {threads} -m {params.mem_mb_per_thread}M -n -o $LOCAL/temp.qsorted.bam -T /local/ $LOCAL/raw.bam
            set +e
            rand_str=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 16)
            set -e
            temp_file_rand="$temp_file.$rand_str"
            cp -av $LOCAL/temp.qsorted.bam "temp/bwa_dc/$temp_file_rand"
            mv "temp/bwa_dc/$temp_file_rand" "temp/bwa_dc/$temp_file"
        fi

        echo `date +"%F %T %Z"` "samblaster processing..."
        samtools view -@ {threads} -h $LOCAL/temp.qsorted.bam | samblaster 2> >(tee $LOCAL/summary.txt >&2) | \
            samtools view -b - > $LOCAL/temp.mdups.bam

        set +e
        rand_str=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 16)
        set -e
        temp_file="temp/bwa_dc/${{output_basename}}.temp.mdups.bam.$rand_str"
        cp -a $LOCAL/temp.mdups.bam $temp_file
        mv $temp_file "temp/bwa_dc/${{output_basename}}.temp.mdups.bam"
        
        echo `date +"%F %T %Z"` "position-sorting..."
        samtools sort -@ {threads} -m {params.mem_mb_per_thread}M -o $LOCAL/output.bam -T /local/ $LOCAL/temp.mdups.bam

        echo `date +"%F %T %Z"` "Indexing..."
        samtools index -@ {threads} $LOCAL/output.bam $LOCAL/output.bai

        echo `date +"%F %T %Z"` "Copying rule results to target directory..."

        set +e
        rand_str=$(tr -dc A-Za-z0-9 < /dev/urandom | head -c 16)
        set -e
        temp_file="temp/bwa_dc/${{output_basename}}.mdups.bam.$rand_str"
        cp -a $LOCAL/output.bam $temp_file
        mv $temp_file {output.bam}

        cp -a $LOCAL/summary.txt {output.summary}
        cp -a $LOCAL/output.bai {output.bai}

        echo `date +"%F %T %Z"` "Cleaning up..."
        rm -f "temp/bwa_dc/${{output_basename}}.temp.qsorted.bam"
        rm -f "temp/bwa_dc/${{output_basename}}.temp.mdups.bam"
        """


rule calc_frag:
    singularity: "docker://zephyre/comp-bio:v0.3.7"
    input:
        sam="temp/filtered_sam/{entry_id}.{assembly}.mdups.sam"
    output:
        "frag/{entry_id,EE[0-9]+}.{assembly,hg[0-9]+}.frag.gz",
        "frag/{entry_id}.{assembly}.frag.gz.tbi"
    threads: 4
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4200,
        time_min=lambda wildcards, input, threads: round(nlogn(os.path.getsize(input.sam) / 3) * 32 / threads + 30),
        input_size=lambda wildcards, input: round(os.path.getsize(input.sam) / 1024**2)
    params:
        mem_mb_per_thread=lambda wildcards, threads, resources: int((resources.mem_mb - 1000) * 0.8 / threads),
        mem_mb_sortbed=lambda wildcards, threads, resources: int((resources.mem_mb - 2000) * 0.4),
        sort_threads=lambda wildcards, threads, resources: resources.cpus - 1,
        label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}",
        partition="RM-shared"
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local

        echo `date +"%F %T %Z"` "Filtering and sorting..."
        cat {input.sam} | \
            samtools sort -@ {params.sort_threads} -n -m {params.mem_mb_per_thread}M -T $LOCAL -o $LOCAL/temp.qsorted.bam -
        
        echo `date +"%F %T %Z"` "Calculating fragments..."
        bamToBed -bedpe -mate1 -i $LOCAL/temp.qsorted.bam | \
            perl -ne 'chomp;@f=split " ";if($f[0] ne $f[3]){{next;}}$s=$f[1];$e=$f[5];if($f[8] eq "-"){{$s=$f[4];$e=$f[2];}}if($e>$s){{print "$f[0]\\t$s\\t$e\\t$f[7]\\t$f[8]\\n";}}' | \
            sort-bed --max-mem {params.mem_mb_sortbed}M --tmpdir $LOCAL/ - | \
            bgzip > $LOCAL/frag.gz

        echo `date +"%F %T %Z"` "Indexing fragments..."
        tabix -0 -p bed $LOCAL/frag.gz

        echo `date +"%F %T %Z"` "Copying results to target directory..."
        cp -a $LOCAL/frag.gz {output[0]}.tmp
        mv {output[0]}.tmp {output[0]}
        cp -a $LOCAL/frag.gz.tbi {output[1]}.tmp
        mv {output[1]}.tmp {output[1]}
        """

rule reads_summary:
    singularity: "docker://zephyre/comp-bio:v0.3.7"
    input:
        bam="bam/{entry_id}.{assembly}.mdups.bam"
    output:
        "reads_summary/{entry_id,EE[0-9]+}.{assembly,hg[0-9]+}.reads_summary.txt"
    threads: 1
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4000
    params:
        label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}",
        partition="RM-shared"
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"

        python scripts/reads_summary.py {input.bam} {output}
        """


rule bam2sam:
    singularity: "docker://zephyre/comp-bio:v0.3.7"
    input:
        bam="bam/{entry_id}.{assembly}.mdups.bam"
    output:
        sam="temp/filtered_sam/{entry_id,EE[0-9]+}.{assembly,hg[0-9]+}.mdups.sam"
    threads: 1
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4200,
        time_min=lambda wildcards, input: math.ceil(os.path.getsize(input.bam) / 1024**3) * 10 + 45,
        input_size=lambda wildcards, input: round(os.path.getsize(input.bam) / 1024**2)
    params:
        label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}",
        partition="RM-shared"
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local

        echo `date +"%F %T %Z"` "Filtering..."
        samtools view -h -f 3 -F 3852 -q 30 {input.bam}> $LOCAL/filtered.sam

        echo `date +"%F %T %Z"` "Copying results to target directory..."
        cp -a $LOCAL/filtered.sam {output.sam}.tmp
        mv {output.sam}.tmp {output.sam} 
        """


rule picard_metrics_insert_size:
    singularity: "docker://zephyre/comp-bio:v0.3.7"
    input:
        sam="temp/filtered_sam/{entry_id}.{assembly}.mdups.sam"
    output:
        insert_size="picard_metrics/insert_size/{entry_id,EE[0-9]+}.{assembly,hg[0-9]+}.insert_size.txt",
        insert_size_pdf="picard_metrics/insert_size/{entry_id}.{assembly}.insert_size.pdf",
        insert_size_log="picard_metrics/insert_size/{entry_id}.{assembly}.insert_size.log"
    threads: 1
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4200,
        time_min=lambda wildcards, input: math.ceil(os.path.getsize(input.sam) / 3 / 1024**3) * 10 + 45,
        input_size=lambda wildcards, input: round(os.path.getsize(input.sam) / 1024**2)
    params:
        label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}",
        partition="RM-shared",
        xmx_mb=lambda wildcards, threads: threads * 4000 - 768,
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local

        echo `date +"%F %T %Z"` "Picard insert size analysis..."
        picard -XX:+UseParallelGC -Dpicard.useLegacyParser=false -Xms512m -Xmx{params.xmx_mb}m -XX:MaxMetaspaceSize=768m \
            CollectInsertSizeMetrics -I {input.sam} -O $LOCAL/output.insert_size.txt -H $LOCAL/output.insert_size.pdf \
            2> >(tee $LOCAL/output.insert_size.log >&2)

        echo `date +"%F %T %Z"` "Copying results to target directory..."
        cp -a $LOCAL/output.insert_size.txt {output.insert_size}.tmp
        mv {output.insert_size}.tmp {output.insert_size}
        cp -a $LOCAL/output.insert_size.pdf {output.insert_size_pdf}.tmp
        mv {output.insert_size_pdf}.tmp {output.insert_size_pdf}
        cp -a $LOCAL/output.insert_size.log {output.insert_size_log}.tmp
        mv {output.insert_size_log}.tmp {output.insert_size_log}
        """

rule picard_metrics_lib_complexity:
    singularity: "docker://zephyre/comp-bio:v0.3.7"
    input:
        sam="temp/filtered_sam/{entry_id}.{assembly}.mdups.sam",
    output:
        lc="picard_metrics/lib_complexity/{entry_id,EE[0-9]+}.{assembly,hg[0-9]+}.lib_complexity.txt",
        lc_log="picard_metrics/lib_complexity/{entry_id}.{assembly}.lib_complexity.log"
    threads: 1
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4200,
        time_min=lambda wildcards, input: round(nlogn(os.path.getsize(input.sam) / 3) * 10 + 45),
        input_size=lambda wildcards, input: round(os.path.getsize(input.sam) / 1024**2)
    params:
        label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}",
        partition="RM-shared",
        xmx_mb=lambda wildcards, threads: threads * 4000 - 768
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local

        echo `date +"%F %T %Z"` "Picard library complexity analysis..."
        picard -XX:+UseParallelGC -Dpicard.useLegacyParser=false -Xms512m -Xmx{params.xmx_mb}m -XX:MaxMetaspaceSize=768m \
            EstimateLibraryComplexity -I {input.sam} -O $LOCAL/output.txt -VALIDATION_STRINGENCY SILENT \
            2> >(tee $LOCAL/output.log >&2)

        echo `date +"%F %T %Z"` "Copying results to target directory..."
        cp -a $LOCAL/output.txt {output.lc}.tmp
        mv {output.lc}.tmp {output.lc}
        cp -a $LOCAL/output.log {output.lc_log}.tmp
        mv {output.lc_log}.tmp {output.lc_log}
        """


rule picard_metrics_gc_bias:
    singularity: "docker://zephyre/comp-bio:v0.3.7"
    input:
        sam="temp/filtered_sam/{entry_id}.{assembly}.mdups.sam",
        ref=lambda wildcards: ref_genome_fa[wildcards.assembly]["path"]
    output:
        gc_bias="picard_metrics/gc_bias/{entry_id,EE[0-9]+}.{assembly,hg[0-9]+}.gc_bias.txt",
        gc_bias_chart="picard_metrics/gc_bias/{entry_id}.{assembly}.gc_bias.pdf",
        gc_bias_summary="picard_metrics/gc_bias/{entry_id}.{assembly}.gc_bias.summary.txt",
        gc_bias_log="picard_metrics/gc_bias/{entry_id}.{assembly}.gc_bias.log"
    threads: 1
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4200,
        time_min=lambda wildcards, input: round(nlogn(os.path.getsize(input.sam) / 3) * 10 + 30),
        input_size=lambda wildcards, input: round(os.path.getsize(input.sam) / 1024**2)
    params:
        xmx_mb=lambda wildcards, threads: threads * 4000 - 768,
        label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}",
        partition="RM-shared"
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local

        echo `date +"%F %T %Z"` "Picard gc bias analysis..."
        picard -XX:+UseParallelGC -Dpicard.useLegacyParser=false -Xms512m -Xmx{params.xmx_mb}m -XX:MaxMetaspaceSize=768m \
            CollectGcBiasMetrics -I {input.sam} -R {input.ref} \
            -O $LOCAL/gc_bias.txt -CHART $LOCAL/gc_bias.pdf -S $LOCAL/gc_bias.summary.txt \
            2> >(tee $LOCAL/gc_bias.log >&2)

        echo `date +"%F %T %Z"` "Copying results to target directory..."
        cp -a $LOCAL/gc_bias.txt {output.gc_bias}.tmp
        mv {output.gc_bias}.tmp {output.gc_bias}
        cp -a $LOCAL/gc_bias.pdf {output.gc_bias_chart}.tmp
        mv {output.gc_bias_chart}.tmp {output.gc_bias_chart}
        cp -a $LOCAL/gc_bias.summary.txt {output.gc_bias_summary}.tmp
        mv {output.gc_bias_summary}.tmp {output.gc_bias_summary}
        cp -a $LOCAL/gc_bias.log {output.gc_bias_log}.tmp
        mv {output.gc_bias_log}.tmp {output.gc_bias_log}
        """


rule calc_coverage:
    singularity: "docker://zephyre/comp-bio:v0.3.7"
    input:
        chrom_sizes=lambda wildcards: chrom_sizes[wildcards.assembly],
        ucsc_hg19=lambda wildcards: chrom_sizes["ucsc_hg19"],
        frag="frag/{entry_id}.{assembly}.frag.gz"
    output:
        # Calculate genome coverage for both MAPQ>=30 and whatever MAPQ
        bw="coverage/{entry_id}.{assembly}.coverage.mapq30.bw"
        # # Output
        # bg=expand("coverage/{{entry_id}}.{{assembly}}.coverage.mapq30.bedGraph.gz{postfix}", postfix=["", ".tbi"])
    threads: 2
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4200,
        time_min=lambda wildcards, input, threads: round(nlogn(os.path.getsize(input.frag)) * 60 / threads + 30),
        input_size=lambda wildcards, input: round(os.path.getsize(input.frag) / 1024**2)
    params:
        mem_mb_sortbed=lambda wildcards, threads, resources: int((resources.mem_mb - 2000) * 0.5),
        label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}",
        partition="RM-shared"
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local

        # MAPQ>=30
        echo `date +"%F %T %Z"` "Calculating the MAPQ>=30 coverage..."
        bgzip -d < {input.frag} | awk '$4>=30' | bedtools genomecov -bg -i - -g {input.chrom_sizes} | bedClip /dev/stdin {input.chrom_sizes} $LOCAL/temp.bedGraph

        if [[ "{wildcards.assembly}" == "hg19" ]]; then
            echo `date +"%F %T %Z"` "Assembly is hg19, will go through contig mapping"
            cat $LOCAL/temp.bedGraph | python scripts/contig_mapping.py scripts/contig_mapping.csv | sort-bed --max-mem {params.mem_mb_sortbed}M --tmpdir $LOCAL/ - > $LOCAL/temp.contig_mapped.bedGraph
            bedGraphToBigWig $LOCAL/temp.contig_mapped.bedGraph {input.ucsc_hg19} $LOCAL/output.bw
        else
            echo `date +"%F %T %Z"` "Assembly is not hg19, contig mapping is not needed"
            bedGraphToBigWig $LOCAL/temp.bedGraph {input.chrom_sizes} $LOCAL/output.bw
        fi

        echo `date +"%F %T %Z"` "Copying results to target directory..."
        cp -a $LOCAL/output.bw {output.bw}.tmp
        mv {output.bw}.tmp {output.bw}
        """

# rule calc_frag_profile:
#     singularity: "docker://zephyre/comp-bio:v0.3.7"
#     input:
#         chrom_sizes=lambda wildcards: chrom_sizes[wildcards.assembly],
#         ucsc_hg19=lambda wildcards: chrom_sizes["ucsc_hg19"],
#         frag="frag/{entry_id}.{assembly}.frag.gz"
#     output:
#         # Calculate genome coverage for both MAPQ>=30 and whatever MAPQ
#         bw="frag_profile/{entry_id}.{assembly}.frag_profile.mapq30.bw",
#         # Output
#         bg=expand("frag_profile/{{entry_id}}.{{assembly}}.frag_profile.mapq30.bedGraph.gz{postfix}", postfix=["", ".tbi"])
#     threads: 2
#     resources:
#         cpus=lambda wildcards, threads: threads,
#         mem_mb=lambda wildcards, threads: threads * 4200,
#         time_min=lambda wildcards, input, threads: round(nlogn(os.path.getsize(input.frag)) * 60 / threads + 30),
#         input_size=lambda wildcards, input: round(os.path.getsize(input.frag) / 1024**2)
#     params:
#         mem_mb_sortbed=lambda wildcards, threads, resources: int((resources.mem_mb - 2000) * 0.5),
#     shell:
#         """
#         echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
#         LOCAL=/local
#         TMPDIR=/local/tmp
#         cp -a {input.ucsc_hg19} $LOCAL/ucsc_hg19.sizes

#         # MAPQ>=30
#         echo `date +"%F %T %Z"` "Calculating the MAPQ>=30 coverage..."
#         bgzip -d < {input.frag} | awk '$4>=30' | bedtools genomecov -bg -i - -g {input.chrom_sizes} | bedClip /dev/stdin {input.chrom_sizes} $LOCAL/temp.bedGraph

#         if [[ "{wildcards.assembly}" == "hg19" ]]; then
#             echo `date +"%F %T %Z"` "Assembly is hg19, will go through contig mapping"
#             cp $LOCAL/temp.bedGraph temp/{wildcards.entry_id}.{wildcards.assembly}.coverage.temp.bedGraph
#             cat $LOCAL/temp.bedGraph | python scripts/contig_mapping.py scripts/contig_mapping.csv | sort-bed --max-mem {params.mem_mb_sortbed}M --tmpdir $LOCAL/ - > $LOCAL/temp.contig_mapped.bedGraph
#             cp $LOCAL/temp.contig_mapped.bedGraph temp/{wildcards.entry_id}.{wildcards.assembly}.coverage.temp.contig_mapped.bedGraph
#             bedGraphToBigWig $LOCAL/temp.contig_mapped.bedGraph {input.ucsc_hg19} $LOCAL/output.bw
#         else
#             echo `date +"%F %T %Z"` "Assembly is not hg19, contig mapping is not needed"
#             bedGraphToBigWig $LOCAL/temp.bedGraph {input.chrom_sizes} $LOCAL/output.bw
#         fi

#         echo `date +"%F %T %Z"` "Compressing coverage bedGraph..."
#         bgzip -@ {threads} $LOCAL/temp.bedGraph

#         echo `date +"%F %T %Z"` "Indexing coverage bedGraph..."
#         tabix -p bed $LOCAL/temp.bedGraph.gz

#         echo `date +"%F %T %Z"` "Copying results to target directory..."
#         parallel -n 2 cp -a ::: "$LOCAL/temp.bedGraph.gz" {output.bg[0]} "$LOCAL/temp.bedGraph.gz.tbi" {output.bg[1]} "$LOCAL/output.bw" {output.bw}

#         echo `date +"%F %T %Z"` "Clearing up..."
#         rm -f temp/{wildcards.entry_id}.{wildcards.assembly}.coverage.*


#         #         # MAPQ>=30
# #         echo "Calculating the MAPQ>=30 fragment profile..."
# #         bgzip -d < {input.frag} | awk '$4>=30' | python calc_frag_profile.mosdepth.py -i - -g {input.chrom_sizes} -w 500 | bedClip /dev/stdin {input.chrom_sizes} temp.bedGraph
# #         bgzip -c temp.bedGraph > {output.bg[0]}
# #         tabix -p bed {output.bg[0]}
# #         if [[ "{wildcards.assembly}" == "hg19" ]]; then
# #             echo "Assembly is hg19, will go through contig mapping"
# #             unlink temp.bedGraph
# #             bgzip -d < {output.bg[0]} | python contig_mapping.py | sort-bed --max-mem {params.mem_mb_sortbed}M - | bedClip /dev/stdin {input.ucsc_hg19} temp.contig_mapped.bedGraph
# #             bedGraphToBigWig temp.contig_mapped.bedGraph {input.ucsc_hg19} {output.bw}
# #             unlink temp.contig_mapped.bedGraph
# #         else
# #             echo "Assembly is not hg19, contig mapping is not needed"
# #             bedGraphToBigWig temp.bedGraph {input.chrom_sizes} {output.bw}
# #             unlink temp.bedGraph
# #         fi

# #         touch {output.log}
#         """

rule calc_wps:
    singularity: "docker://zephyre/comp-bio:v0.3.7"
    input:
        chrom_sizes=lambda wildcards: chrom_sizes[wildcards.assembly],
        ucsc_hg19=lambda wildcards: chrom_sizes["ucsc_hg19"],
        frag="frag/{entry_id}.{assembly}.frag.gz"
    output:
        bw="wps/{entry_id}.{assembly}.{wps_type}.mapq30.bw",
    threads: 1
    resources:
        cpus=lambda wildcards, threads: threads,
        mem_mb=lambda wildcards, threads: threads * 4200,
        time_min=lambda wildcards, input, threads: round(nlogn(os.path.getsize(input.frag)) * 90 / threads + 60),
        input_size=lambda wildcards, input: round(os.path.getsize(input.frag) / 1024**2)
    params:
        mem_mb_sortbed=lambda wildcards, threads, resources: int((resources.mem_mb - 2000) * 0.5),
        label=lambda wildcards: f"{wildcards.entry_id}.{wildcards.assembly}",
        partition="RM-shared"
    shell:
        """
        echo `date +"%F %T %Z"` "Job started: entry_id: {wildcards.entry_id}"
        LOCAL=/local
        TMPDIR=/local/tmp
        cp -a {input.ucsc_hg19} $LOCAL/ucsc_hg19.sizes

        # MAPQ>=30
        echo `date +"%F %T %Z"` "Calculating the MAPQ>=30 WPS..."
        if [[ "{wildcards.wps_type}" == "l_wps" ]]; then
            bgzip -d < {input.frag} | awk '$4>=30 && $3-$2>120' | python scripts/calc_wps.py -i - -g {input.chrom_sizes} -w 60 | bedClip /dev/stdin {input.chrom_sizes} $LOCAL/temp.bedGraph
        elif [[ "{wildcards.wps_type}" == "s_wps" ]]; then
            echo "Invalid WPS type: {wildcards.wps_type}"
            exit 1
        else
            echo "Invalid WPS type: {wildcards.wps_type}"
            exit 1
        fi

        if [[ "{wildcards.assembly}" == "hg19" ]]; then
            echo `date +"%F %T %Z"` "Assembly is hg19, will go through contig mapping"
            cat $LOCAL/temp.bedGraph | python scripts/contig_mapping.py scripts/contig_mapping.csv | sort-bed --max-mem {params.mem_mb_sortbed}M --tmpdir $LOCAL/ - > $LOCAL/temp.contig_mapped.bedGraph
            echo `date +"%F %T %Z"` "Generating the bigWig track..."
            bedGraphToBigWig $LOCAL/temp.contig_mapped.bedGraph {input.ucsc_hg19} $LOCAL/output.bw
        else
            echo `date +"%F %T %Z"` "Assembly is not hg19, contig mapping is not needed"
            echo `date +"%F %T %Z"` "Generating the bigWig track..."
            bedGraphToBigWig $LOCAL/temp.bedGraph {input.chrom_sizes} $LOCAL/output.bw
        fi

        echo `date +"%F %T %Z"` "Copying results to target directory..."
        cp -a $LOCAL/output.bw {output.bw}.tmp
        mv {output.bw}.tmp {output.bw}
        """


# rule samtools_stats:
#     input:
#         bam="entries/{entry_id}/{assembly}/{entry_id}.{assembly}.mdups.bam"
#     output:
#         "entries/{entry_id}/{assembly}/{entry_id}.{assembly}.stats.txt",
#         "entries/{entry_id}/{assembly}/{entry_id}.{assembly}.filtered.stats.txt",
#     threads: 2
#     shell:
#         """
#         samtools stats -@ {threads} {input.bam} > {output[0]}
#         samtools stats -@ {threads} -f 3 -F 3852 {input.bam} > {output[1]}
#         """
