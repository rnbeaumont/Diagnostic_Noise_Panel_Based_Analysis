#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
my $list=undef;
my $nquant=5;
my $outname=undef;
my $ukb_dir="./"
GetOptions(
	"list=s"     => \$list,
	"out=s"      => \$outname,
	"quant=s"    => \$nquant,
	"ukb_dir=s"  => \$ukb_dir,
);
if($list && $nquant && $outname){
	# first read in the list of all SNPs (do I want quantiles per gene or across all?)
	my @revel=();
	my %snp=();
	for(my $i=1;$i<23;$i++){
		open(my $in,"results/".$i."_G2P_".$list.".txt") or die $!;
		while(<$in>){
			chomp;
			my @F=split("\t");
			my @G=split(";",$F[6]);
			for(my $j=0;$j<@G;$j++){
				my @H=split(":",$G[$j]);
				$snp{$H[0].":".$H[1]}=1;
			}
		}
		close($in);
	}
	# now read in their revel scores from the VEP annotation Andy did
	open(my $out,">",$outname);
	for(my $i=1;$i<23;$i++){
		open(my $in,"zcat ${ukb_dir}/UKBexomeOQFE_chr".$i."_vep.gz | ") or die $!;
		my $line=undef;
		do{
			$line=<$in>;
		}while($line !~ /^#Uploaded_variation/);
		while(<$in>){
			chomp;
			my @F=split(' ');
			if(exists $snp{$F[1]}){
				my @G=split(";",$F[13]);
				for(my $j=0;$j<@G;$j++){
					my @H=split("=",$G[$j]);
					if($H[0] eq "REVEL"){
						push @revel, $H[1];
						print $out $F[1]."\t".$F[2]."\t".$H[1]."\n";
					}
				}
			}
		}
		close($in);
	}
	close($out);
	# now sort the vector
	@revel=sort @revel;
	# then print quantiles
	my $n=@revel;
	my $cut=$n/$nquant;
	print "Quantile cutoffs:\n";
	for(my $i=0;$i<$nquant-1;$i++){
		print $revel[int($cut*($i+1))]."\n";
	}
	print $revel[@revel]."\n";
}else{
	print "\n\tUsage:\n";
	print "\t--list         Gene list name\n";
	print "\t--out          Output filename\n";
	print "\t--nquant       (optional) Number of quantiles to calculate. Default 5\n\n";
}
