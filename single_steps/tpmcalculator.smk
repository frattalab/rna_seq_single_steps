import os
import subprocess

include: "helpers.py"

configfile: "config/tpmcalculator.yaml"

import os

project_folder = config["project_top_level"]
end_type = config["end_type"]
suffix = config["suffix"]
bam_folder = project_folder + config['bam_subfolder']
tpm_output_folder = project_folder + config["tpmcalculator_output_folder"]
REFERENCE_ANNOTATION = config["reference_gtf"]
SAMPLE_NAMES, = glob_wildcards(bam_folder + "{sample}" + suffix + ".bam")

print(SAMPLE_NAMES)

#this function uses the text file located in the config folder "star_genomes_species.csv" and
#the config file species parameter to
#give the correct genome for the species
os.system("mkdir -p {0}".format("tpm_output_folder"))
os.system("cd {0}".format("tpm_output_folder"))


rule tpmcounts:
    input:
        expand(tpm_output_folder + "{name}" + suffix + "_genes.out", name = SAMPLE_NAMES)


rule tpmcalculator_path:
    input:
        aligned_bam = bam_folder + "{name}" + suffix + ".bam",
        aligned_bai = bam_folder + "{name}" + suffix + ".bam.bai"
    output:
        tpm_output_folder + "{name}" + suffix + "_genes.out"
    params:
        ref_anno = REFERENCE_ANNOTATION,
        output_string = "-p -e -a" if config["end_type"] == "pe" else "-e -a"
    shell:
        """
        cd {tpm_output_folder}
        pwd
        {config[tpmcalculator_path]} -g {params.ref_anno} -b {input.aligned_bam} {params.output_string}
        """
