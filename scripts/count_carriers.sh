#!/bin/bash

PREFIX=$1
chr=$2
dir=$3
./scripts/count_ukbb_carriers.pl --vep ./results/${chr}_G2P_${PREFIX}.txt --plink ${dir}/ukb23155_c${chr}_b0_v1 --out ./results/chr_${chr}_variants_${PREFIX} --freq ./gnomad_files/gnomad_chr${chr}_freq --freq_ukb ./results/${chr}_freq.afreq
