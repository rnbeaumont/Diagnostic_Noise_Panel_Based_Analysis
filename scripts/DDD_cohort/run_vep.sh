#!/bin/bash

directory=$1
mkdir results/
#first run list_variants.sh to get list of variants
./list_variants.pl

# then run VEP on the lsit of variants
vep -i DDD_variant_list -o results/DDD_variants_annotated_orig --cache --assembly GRCh37 --plugin G2P,file="../g2p_files/DDG2P.csv,log_dir=results/DDD_variants_annotated_G2P_report_orig,html_report=results/DDD_variants_annotated_G2P_orig.html,txt_report=results/DDD_variants_annotated_G2P_orig.txt"

# count the variants
./count_ddd_carriers.pl --vep ./results/DDD_variants_annotated_G2P.txt --plink ${directory} --out results/DDD_all_chr_variant_counts

# Summarise genes and individuals
./sum_variants.pl --infile results/DDD_all_chr_variant_counts --trios ${directory}/trios.txt --outfile results/DDD_indivudals_variants
./sum_variants.pl --infile results/DDD_all_chr_variant_counts --trios DDD_parents --outfile results/DDD_indivudals_variants_parents
./sum_variants_missense_lof.pl --infile results/DDD_all_chr_variant_counts_missense_lof --trios /slade/projects/DDD/PTB_data/trios.txt --outfile results/DDD_indivudals_variants_missense_lof

# then examine biallelic genes for cis vs trans
./count_ddd_carriers_biallelic_cis_trans.pl --vep results/DDD_variants_annotated_G2P.txt --vep_rdif results/DDD_RD_IF_variants_annotated_G2P.txt --plink $directory --out results/DDD_all_chr_variant_counts_biallelic_cis_trans --gene_counts results/DDD_all_chr_variant_counts.counts
./sum_genes_cis_trans.pl --in results/DDD_all_chr_variant_counts_biallelic_cis_trans --out results/DDD_genes_counts_cis_trans --individuals /slade/projects/DDD/PTB_data/trios.txt
