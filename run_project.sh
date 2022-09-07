#!/bin/bash
# script to run entire project
# The project will require access to UKBB and DDD data and VEP with the G2P-VEP plugin installed. Cache files for GRCh37 and GRCh38 are also required

ukb_dir="./"	# directory containing UKBB data
g2p_dir="./"	# directory containing G2P panel files
ddd_dir="./DDD_cohort/"	# directory containing DDD data

mkdir results/

# reformat UKBB data to make dummy VCF file
for i in {1..22} X Y
do
	awk 'NR==1{print;next} NR==2{print $0"\tQUAL\tFILTER\tINFO\tFORMAT\tID_1\tID_2\tID_3";next} {print $0"\t42.5\tPASS\t.\tGT:GQ:MIN_DP\t0/0:30:30\t0/1:30:30\t1/1:30:30"}' ${ukb_dir}UKBexomeOQFE_chr${i}.vcf | bgzip -c > UKBexomeOQFE_chr${i}_dummy.vcf.gz
done

# calculate UKB allele frequency
generate_af () {
        plink2 --bfile ${ukb_dir}/ukb23155_c${1}_b0_v1 --freq --out results/${1}_freq
}
export -f generate_af
echo {1..22} | xargs -n 1 -P 12 bash -c 'generate_af "$@"' _

# count variants using wrapper scripts
for list in DD Eye Cancer Skin Cardiac
do
	./run_vep.sh $list
	for chr in {1..22}
	do
		./count_carriers.sh $list $chr $ukb_dir
	done
done

# then loop over the gene lists and count the carriers for each
for PREFIX in DD Eye Cancer Skin Cardiac
do
	# calculate REVEL quintiles
	./calculate_revel_quantiles_${PREFIX}.sh

	# then calculate total variants by summing over all chromosomes and all individuals
	./sum_variants.pl --infile ${PREFIX} --outfile results/${PREFIX}_individual_counts_chrs_collated	# by individual
	./sum_genes.pl --in ${PREFIX} --out results/${PREFIX}_gene_variant_counts --glist g2p_files/${PREFIX}G2P.csv	# by gene
	./sum_variants_missense_lof.pl --infile ${PREFIX}_missense_lof --outfile results/${PREFIX}_individual_counts_chrs_collated_missense_lof --vep_rep ${PREFIX}	# missense and lof split

	./sum_genes_all_transcripts.pl --in ${PREFIX}_all_transcripts --out results/${PREFIX}_all_transcripts_gene_variant_counts --glist g2p_files/${PREFIX}G2P.csv
	./sum_genes_missense_lof.pl --in ${PREFIX}_missense_lof --out results/${PREFIX}_missense_lof_gene_variant_counts --glist g2p_files/${PREFIX}G2P.csv
	# clinvar
	./sum_genes_clinvar.pl --in ${PREFIX}_clinvar --out results/${PREFIX}_clinvar_gene_variant_counts --glist g2p_files/${PREFIX}G2P.csv
	./sum_variants_clinvar.pl --infile ${PREFIX}_clinvar --outfile results/${PREFIX}_individual_counts_chrs_collated_clinvar --vep_rep ${PREFIX}
	# revel
	./sum_genes_revel.pl --in ${PREFIX}_revel --out results/${PREFIX}_revel_gene_variants_counts --glist g2p_files/${PREFIX}G2P.csv
	# clinvar and revel thresholds
	./sum_variants_clinvar.pl --infile ${PREFIX}_clinvar_revel --outfile results/${PREFIX}_individual_counts_chrs_collated_clinvar_revel --vep_rep ${PREFIX}
done

PREFIX=Cancer
./sum_genes.pl --in ${PREFIX} --out results/${PREFIX}_cases_gene_variant_counts --individuals cancer_cases --glist g2p_files/${PREFIX}G2P.csv
./sum_genes.pl --in ${PREFIX} --out results/${PREFIX}_controls_gene_variant_counts --individuals cancer_controls --glist g2p_files/${PREFIX}G2P.csv

./DDD_cohort/run_vep.sh $ddd_dir

Rscript render.R
