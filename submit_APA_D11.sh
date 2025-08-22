#!/bin/bash

cd /scratch/nautilus/projects/CRCI2NA_DATA/CHILD/RNAseq_APA_D11
module load guix/latest/
guix time-machine -C .guix/snakemake_7.7.0/channels.scm -- shell -m .guix/snakemake_7.7.0/manifest.scm -- \
snakemake --config projf="project.json" conff="config.json" --jobs 60 --rerun-incomplete --resources parallel_star=3 --resources mem_mb=80000
