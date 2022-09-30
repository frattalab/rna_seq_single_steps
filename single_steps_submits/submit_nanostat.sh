#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=20:00:00,tmem=1.9G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y

if [ "$1" != "" ]; then
    RUN_NAME=$1
else
    RUN_NAME=$""
fi

FOLDER=$(date +"%Y%m%d%H%M")
WRITEFOLDER=../submissions/$FOLDER
mkdir -p $WRITEFOLDER
cp single_steps/nanostat.smk $WRITEFOLDER/nanostat.smk

snakemake -s single_steps/nanostat.smk \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster/nanostat.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o $FOLDER {cluster.submission_string}" \
-j 40 \
--nolock \
--rerun-incomplete \
--latency-wait 100 \
--use-singularity \
--singularity-args "-B /SAN/vyplab/"
