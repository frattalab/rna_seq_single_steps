configfile: "config/samtools_stats_config.yaml"

import os
import sys

in_bam_dir = config["input_dir"]
bam_suffix = config["bam_suffix"]
out_dir = config["output_dir"]

summary_out = os.path.join(out_dir, config["summary_out_file"])
log_dir = os.path.join(out_dir, config["log_dir"])

SAMPLES = [f.replace(bam_suffix, "") for f in os.listdir(in_bam_dir) if f.endswith(bam_suffix)]

# Check that each input sample has a BAM index (.bai)
for s in SAMPLES:
    assert os.path.exists(os.path.join(in_bam_dir, s + bam_suffix + ".bai")), f".bai index file does not exist at same location as input BAM file for sample {s}"

assert isinstance(config["remove_pe_overlaps"], bool), f"'remove_pe_overlaps' must be True/False boolean, {config['remove_pe_overlaps']} (type {type(config['remove_pe_overlaps'])}) was provided"


sys.stderr.write(f"Basenames for input BAM files - {', '.join(SAMPLES)}\n")

if not os.path.exists(out_dir):
    os.system(f"mkdir -p {out_dir}")


wildcard_constraints:
    sample = "|".join(SAMPLES)


rule all:
    input:
        # summary_out,
        expand(os.path.join(out_dir, "{sample}.samtools_stats.txt"), sample=SAMPLES)


rule samtools_stats:
    input:
        os.path.join(in_bam_dir, "{sample}" + bam_suffix)

    output:
        os.path.join(out_dir, "{sample}.samtools_stats.txt")

    params:
        extra_threads = config["stats_extra_threads"],
        rm_overlaps = "--remove-overlaps" if config["remove_pe_overlaps"] else "",

    threads:
        config["stats_extra_threads"] + 1

    log:
        os.path.join(log_dir, "{sample}.samtools_stats.log")

    conda:
        "../envs/single_steps.yaml"

    shell:
        """
        samtools stats \
        --threads {params.extra_threads} \
        {params.rm_overlaps} \
        {input} > {output} 2> {log}
        """
