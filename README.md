# Snakemake workflows for FinaleDB

Snakemake supports running workflows in a number of computing environments, such as local, Kubernetes, SLURM clusters, to name a few. Our workflow are tailored for both Kubernetes and SLURM.

## Work with SLURM batch system

We suggest invoking snakemake workflows with `--profile` (https://snakemake.readthedocs.io/en/stable/executing/cli.html?highlight=profile#profiles). For more information, please refer to [this excellent blog post](https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/). Below is a profile configiuration for example:


```yaml
jobs: 1000
cluster: "sbatch -t 4320 -p {params.partition} --mem={resources.mem_mb} --time-min {resources.time_min} --ntasks-per-node {resources.cpus} --job-name snakemake.{rule}.{params.label} -o slurm_logs/slurm.jobid_%j.{rule}.{params.label}.log"
default-resources: [cpus=1, time_min=1440, mem_mb=4096]
singularity-prefix: "direcotry for keeping singularity images"
singularity-args: "--bind $LOCAL:/local"
directory: "your working directory"
keep-going: true
```

In the above example, snakemake scheduler will submit batch jobs that:
* Request running time capped by 3 days (4320 minutes)
* Request memory of {resources.mem_mb} MB specified by each rule
* Request minimum running time of {resources.time_min} minutes specified by each rule
* Request {resources.cpus} cores
* By default, each rule needs 1 core, 4GB memory, and at 1440 minutes. This can be overridden by the {resources} directive in the rule specification

Note: `directory` specifies the working directory and `-o slurm_logs/...` tells SLURM where to put log files. The log file directory is relative to the working directory. You have to make sure `$directory/slurm_logs` exists prior to run the pipeline.

An example of invoking snakemake:

```bash
snakemake --use-singularity --profile ~/.config/snakemake/finaledb/ -s workflow.slurm.smk --reason -np bam/EE00001.mdups.bam
```

## Work with SLURM interactive nodes

It's fairly similar to working with SLURM batch. You just need to modify the profile configuration by removing the `cluster` directive:

```yaml
jobs: 1000
default-resources: [cpus=1, time_min=1440, mem_mb=4096]
singularity-prefix: "direcotry for keeping singularity images"
singularity-args: "--bind $LOCAL:/local"
directory: "your working directory"
keep-going: true
```

When invoking snakemake from command line, do remember to correctly specify the number of cores. This can be down by `--config FULL_CORES={# of cores}`.

## Run locally

Similar to working with SLURM interactive nodes, except for a fundamental difference: SLURM nodes have node-local storage at $LOCAL (/local). Some of the rules use this directory to keep intermediate output files. However, when you run snakemake locally, there is a good chance `/local` doesn't exist. In this case, you need to modify the profile configuration with a new singularity filesystem mount rule:

```
singularity-args: "--bind $LOCAL:{absolute path for intermediate output}"
```