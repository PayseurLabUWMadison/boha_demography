#!/bin/sh
VCF=$1
KEEP=$2
OUT=$3
# Evalutes singletons for *specified* population, excluding sites with missing data
vcftools --gzvcf ${VCF} --keep ${KEEP} --max-missing 1.0 --singletons --out ${OUT}