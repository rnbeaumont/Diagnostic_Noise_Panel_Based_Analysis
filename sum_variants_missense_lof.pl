#!/usr/bin/perl
# script to collate variant counts into a single sum
use warnings;
use strict;
use Getopt::Long;
my $gene_list=undef;
my $outfile=undef;
my $vep_rep=undef;
GetOptions(
	"infile=s"    => \$gene_list,
	"vep_rep=s"   => \$vep_rep,
	"outfile=s"   => \$outfile,
);
my %count_variants=();
my %count=();
my %count_genes=();
my %count_homs=();
my %variants_homs=();
my %variant_info=();
my %lof=();
my %mis=();
my %low=();
my %mod=();
my %unk=();

if($gene_list && $outfile && $vep_rep){
my $in=undef;
for(my $i=1;$i<23;$i++){	# loop over the autosomes
	open($in,"results/chr_${i}_variants_${gene_list}.counts") or die $!;
	my $line=<$in>;	# skip the header
	chomp $line;
	my @H=split(' ',$line);
	while(<$in>){
		my %ghash=();	# add all genes to a hash so I can count them properly
		chomp;
		my @F=split(' ');
		$count_variants{$F[0]}+=$F[1];
		for(my $i=2;$i<@F;$i++){
			$ghash{$H[$i]}+=$F[$i];
			#if($F[$i]>0){
			#	$count{$F[0]}+=$F[$i];
			#	$count_genes{$F[0]}++;
			#}else{
			#	$count{$F[0]}+=0;
			#	$count_genes{$F[0]}+=0;
			#}
                }
		if(!exists $lof{$F[0]}){
			# initialise hashes
			$lof{$F[0]}=0;
			$mis{$F[0]}=0;
			$low{$F[0]}=0;
			$mod{$F[0]}=0;
			$unk{$F[0]}=0;
			$count_variants{$F[0]}=0;
			$count{$F[0]}=0;
			$count_genes{$F[0]}=0;
		}
		for(my $i=2;$i<@H;$i+=4){	# now loop over all the genes, and we can pull in the mis and lof counts together
			my @G=split("_",$H[$i]);
			my $name=$G[0]."_".$G[1]."_".$G[2];
			if($G[2] eq "monoallelic"){
				if($ghash{$name."_lof"}>0){
					$lof{$F[0]}+=$ghash{$name."_lof"};
				}elsif($ghash{$name."_missense"}>0){
					$mis{$F[0]}+=$ghash{$name."_missense"};
				}elsif($ghash{$name."_low"}>0){
					$low{$F[0]}+=$ghash{$name."_low"};
				}elsif($ghash{$name."_modifier"}>0){
					$mod{$F[0]}+=$ghash{$name."_modifier"};
				}
				$count{$F[0]}+=$ghash{$name."_lof"}+$ghash{$name."_missense"}+$ghash{$name."_low"}+$ghash{$name."_modifier"};
				if(($ghash{$name."_lof"}+$ghash{$name."_missense"}+$ghash{$name."_low"}+$ghash{$name."_modifier"})>0){
					$count_genes{$F[0]}++;
				}
			}elsif($G[2] eq "biallelic"){
				if($ghash{$name."_lof"}>1){
					$lof{$F[0]}+=$ghash{$name."_lof"};
				}elsif($ghash{$name."_missense"}>1){
					$mis{$F[0]}+=$ghash{$name."_missense"}+$ghash{$name."_lof"};
				}elsif($ghash{$name."_low"}>1){
					$low{$F[0]}+=$ghash{$name."_low"}+$ghash{$name."_missense"}+$ghash{$name."_lof"};
				}elsif($ghash{$name."_modifier"}>1){
					$mod{$F[0]}+=$ghash{$name."_modifier"}+$ghash{$name."_low"}+$ghash{$name."_missense"}+$ghash{$name."_lof"};
				}
				$count{$F[0]}+=$ghash{$name."_lof"}+$ghash{$name."_missense"}+$ghash{$name."_low"}+$ghash{$name."_modifier"};
				if(($ghash{$name."_lof"}+$ghash{$name."_missense"}+$ghash{$name."_low"}+$ghash{$name."_modifier"})>1){
					$count_genes{$F[0]}++;
				}
			}else{
				$count{$F[0]}+=0;
				$count_genes{$F[0]}+=0;
			}
		}
	}
	close($in);
	open($in,"results/chr_${i}_variants_${gene_list}.variant_info") or die $!;
	<$in>;
	while(<$in>){
		chomp;
		my @F=split(' ',$_,2);
		if(!exists $variant_info{$F[0]}){
			$variant_info{$F[0]}=$F[1];
		}else{
			$variant_info{$F[0]}.="\t".$F[1];
		}
		# now look through the variants to identify homozygous calls
		my @G=split(",",$F[1]);
		for(my $i=0;$i<@G;$i++){
			my @H=split("_",$G[$i]);
			if($H[1]>1){
				$count_homs{$F[0]}++;
				$variants_homs{$F[0]}.="\t".$H[0];
			}
		}
	}
	close($in);
}
open(my $out,">",$outfile.".counts") or die $!;
print $out "ID\tcount_variants_unique\tcount_variants\tcount_genes\tcount_homs\tlof\tmissense\tlow\tmodifier\tunknown\n";
foreach my $key (keys %count){
	print $out join("\t",$key,$count_variants{$key},$count{$key},$count_genes{$key});
	if(exists $count_homs{$key}){
		print $out "\t".$count_homs{$key};
	}else{
		print $out "\t0";
	}
	print $out join("\t","\t".$lof{$key},$mis{$key},$low{$key},$mod{$key},$unk{$key}."\n");
}
close($out);
# now find the genes associated with the homozygote variants
# loop through the VEP results looking for the genes associated with the variants
my %snp_to_gene=();
for(my $i=1;$i<23;$i++){
	my ($snp_hash,$genes,$genes_hash)=&ReadVEPReport("results/".$i."_G2P_".$vep_rep.".txt");
	foreach my $key (%$snp_hash){
		$snp_to_gene{$key}=$snp_hash->{$key};
	}
}
# the output all variants
open($out,">",$outfile.".variant_info") or die $!;
print $out "ID\tvariants\n";
foreach my $key (keys %variant_info){
	my @F=split(' ',$variant_info{$key});
	print $out $key;
	for(my $i=0;$i<@F;$i++){
		my @G=split(",",$F[$i]);
		for(my $j=0;$j<@G;$j++){
			print $out "\t".$G[$j];
			my @H=split("_",$G[$j]);
			if(exists $snp_to_gene{$H[0]}){
				print $out ";".$snp_to_gene{$H[0]};
			}
		}
	}
	print $out "\n";
}
close($out);
# then output the list of variants with gene info attached
open($out,">",$outfile.".variant_info_homozygotes") or die $!;
print $out "ID\tvariants\n";
foreach my $key (keys %variants_homs){
	my @F=split(' ',$variants_homs{$key});
	print $out $key;
	for(my $i=0;$i<@F;$i++){
		print $out "\t".$F[$i];
		if(exists $snp_to_gene{$F[$i]}){
			print $out ";".$snp_to_gene{$F[$i]};
		}
	}
	print $out "\n";
}
close($out);
}else{
	print "\n\tUsage:\n";
	print "\t--infile\tgene list name\n";
	print "\t--vep_rep\tgene list vep report suffix\n";
	print "\t--outfile\tprefix of output files\n\n";
}

sub ReadVEPReport{      # Read in the VEP report and save the results to a hash to return
	    my $vep          = $_[0];
	    my $file_line    = undef;
	    my @data         = ();
	    my @chr          = ();
	    my %snp_hash     = ();
	    my @genes        = ();
	    my %genes_hash   = ();

	    open(my $snp_file, $vep) or die "Unable to open file $vep\n";       # open the file
	    while($file_line=<$snp_file>){      # loop through the lines
	            chomp $file_line;
	            my @F=split("\t",$file_line);   # split out the fields
	            if($F[3] eq "is_canonical"){    # check whether it's labelled as in the canonical transcript
			    my @J=split('=',$F[5]);     # split REQ field
			    my $biallelic = $J[1];
			    $genes_hash{$F[1]."_".$F[2]}=$biallelic;
			    my @H=split(";",$F[6]);     # split up the variants in this gene and loop through them
			            for(my $i=0;$i<@H;$i++){
			                    my @G=split(":",$H[$i]);        # split out the fields in the variant annotation field
			                    if($G[0] eq "X"){       # recode X and Y chr codes so they match VEP output
			                        $G[0]="23";
			                    }elsif($G[0] eq "Y"){
					        $G[0]="24";
					    }
					    # now save the variant info
					    # First we need to check if it's an indel. If so we'll need to fix the positions so they match the plink file
					    if($G[4] eq "-"){       # is's a deletion
			                        my $pos=$G[1]-1;
			                        my @I=split('',$G[3]);      # work out how many bases are deleted
				                my $ndel=@I;
				                if(!exists $snp_hash{$G[0].":".$pos.":D:".$ndel}){
				                    $snp_hash{$G[0].":".$pos.":D:".$ndel}=$F[1]."_".$F[2];
						}
					    }elsif($G[3] eq "-"){
					        my $pos=$G[1]-1;
					        my @I=split('',$G[3]);      # work out how many bases are inserted
					        my $nins=@I;
					        if(!exists $snp_hash{$G[0].":".$pos.":I:".$nins}){
					            $snp_hash{$G[0].":".$pos.":I:".$nins}=$F[1]."_".$F[2];
						    $nins--;
						    $snp_hash{$G[0].":".$G[1].":I:".$nins}=$F[1]."_".$F[2];
				                }
				            }else{  # else it's a SNP, so we can just save it
					        if(!exists $snp_hash{$G[0].":".$G[1].":".$G[3].":".$G[4]}){
					            $snp_hash{$G[0].":".$G[1].":".$G[3].":".$G[4]}=$F[1]."_".$F[2]; # so save whether it's het or hom to the hash
				                }
				            }
			            }
	            }
            }
	    close $snp_file;
	    foreach my $key (keys %genes_hash){
	        push @genes, $key;
            }
	    return(\%snp_hash,\@genes,\%genes_hash);
}
