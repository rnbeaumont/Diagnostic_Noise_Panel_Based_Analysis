#!/bin/bash

module load VEP/104.3-foss-2019a-Perl-5.28.1

SUFFIX=$1

# run VEP
run_vep () {
  vep -i results/UKBexomeOQFE_chr${1}_dummy.vcf.gz -o results/UKBexomeOQFE_chr${1}_G2P_${3} --cache --assembly GRCh38 --plugin G2P,file="${2},log_dir=results/${1}_G2P_report_${3},html_report=results/${1}_G2P_${3}.html,txt_report=results/${1}_G2P_${3}_cahce.txt"
  gzip  results/UKBexomeOQFE_chr${1}_G2P_${SUFFIX}
}

com=""
for i in {1..22};do com=$com$i" g2p_files/"$SUFFIX"G2P.csv "$SUFFIX" ";done
export -f run_vep
echo $com | xargs -n 3 -P 3 bash -c 'run_vep "$@"' _
