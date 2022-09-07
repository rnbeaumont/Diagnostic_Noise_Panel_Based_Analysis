#!/usr/bin/perl
# script to collate variant counts into a single sum by genes rather than individuals
use warnings;
use strict;
use Getopt::Long;
my %count=();
my $infile=undef;
my $outfile=undef;
my $glist=undef;
my $to_include="/slade/projects/Research_Project-MRC158833/UKBiobank/Exome_Data/200k/ukb23155_c22_b0_v1.fam";
GetOptions(
	"in=s"    => \$infile,
	"out=s"   => \$outfile,
	"individuals=s" => \$to_include,
        "glist=s" => \$glist,
);
if($infile && $outfile && $glist){
	&MainSub($infile,$outfile,$to_include,$glist);
}else{
	print "\n\tUsage:\n";
	print "\tscript to sum variants per gene from VEP annotation output\n";
	print "\t--in\tgene panel name\n";
	print "\t--out\toutput filename\n";
	print "\t--individuals\tlist of individuals to include (optional)\n";
        print "\t--glist\tgene list\n\n";
}

sub MainSub{
	my $infile      = $_[0];
	my $outfile     = $_[1];
	my $individuals	= $_[2];
        my $glist       = $_[3];
	my $g_pos       = "ensembl_gene_transcript_lengths";
	my $g_pos_total = "ensembl_gene_positions";
	# first read in the gene lengths
	# first the transcript length
	open(my $in,$g_pos) or die $!;
	<$in>;
	my %g_length=();
	my %gnames=();
	while(<$in>){
		chomp;
		my @F=split("\t");
		my $length=$F[4];
		if($F[5]){
			$length=$F[5];
		}
		$g_length{$F[0]."_".$F[2]}=$length;
		if($F[10]){
			if(exists $gnames{$F[10]}){
				if($length>$g_length{$gnames{$F[10]}}){
					$gnames{$F[10]}=$F[0]."_".$F[2];
				}
			}else{
				$gnames{$F[10]}=$F[0]."_".$F[2];
			}
		}
	}
	close($in);
	# then total length length
	open($in,$g_pos_total) or die $!;
	<$in>;
	my %g_length_tot=();
	while(<$in>){
		chomp;
		my @F=split(' ');
		my $length=$F[5]-$F[4];
		$g_length_tot{$F[0]."_".$F[2]}=$length;
	}
	close($in);
	open($in,$glist) or die $!;
	<$in>;
	my %panel_genes=();
	while(<$in>){
		chomp;
		my @F=split("\"");      # need to be carful splitting out the line as there are commas inside some of the fields!
		my @G=split(",",$F[0]);
		for(my $i=1;$i<@F;$i++){
			if($i%2==0){
				my @H=split(",",$F[$i]);
				for(my $j=1;$j<@H;$j++){
					push @G,$H[$j];
				}
			}else{
				push @G, $F[$i];
			}
		}
		if($G[4] ne "possible" && $G[4] ne "both RD and IF"){
			if(exists $panel_genes{$G[0]}){
				$panel_genes{$G[0]}.=",".$G[5];
			}else{
				$panel_genes{$G[0]}=$G[5];
			}
		}
	}
	close($in);
	# and read in the list of individuals to include
	open($in,$individuals) or die $!;
	my %to_include=();
	while(<$in>){
		chomp;
		my @F=split(' ');
		$to_include{$F[0]}=1;
	}
	close($in);
	# then open the output file and write the header
	open(my $out,">",$outfile);
	open(my $out_het,">",$outfile."_by_type");
	print $out "gene_name\tgene_id\ttranscript_id\tinheritance\tdeleteriousness\tvariants\tindividuals\tlength\ttotal_length\n";
	print $out_het "gene_name\tgene_id\ttranscript_id\tinheritance\tmis_mis\tmis_lof\tlof_lof\tlength\ttotal_length\n";
	# and loop through the chromosomes, reading the variant info for each
	my %genes_found=();
	for(my $i=1;$i<23;$i++){
		ReadFile($infile,$out,$i,\%g_length,\%g_length_tot,\%to_include,$out_het,\%genes_found,\%gnames);
	}
	# now print out 0 for all genes not found
	foreach my $key (keys %panel_genes){
		my @F=split(",",$panel_genes{$key});
		for(my $i=0;$i<@F;$i++){
			if(!exists $genes_found{$key."_".$F[$i]} && ($F[$i] eq "monoallelic" || $F[$i] eq "biallelic")){
				print $out join("\t",$key,"NA","NA",$F[$i],"low","0","0",$g_length{$gnames{$key}},$g_length_tot{$gnames{$key}})."\n"; # need gene lengths here - this requires ensembl ID and transcript ID
				print $out join("\t",$key,"NA","NA",$F[$i],"modifier","0","0",$g_length{$gnames{$key}},$g_length_tot{$gnames{$key}})."\n"; # need gene lengths here - this requires ensembl ID and transcript ID
				print $out join("\t",$key,"NA","NA",$F[$i],"missense","0","0",$g_length{$gnames{$key}},$g_length_tot{$gnames{$key}})."\n"; # need gene lengths here - this requires ensembl ID and transcript ID
				print $out join("\t",$key,"NA","NA",$F[$i],"lof","0","0",$g_length{$gnames{$key}},$g_length_tot{$gnames{$key}})."\n"; # need gene lengths here - this requires ensembl ID and transcript ID
				print $out_het join("\t",$key,"NA","NA",$F[$i],"0","0","0",$g_length{$gnames{$key}},$g_length_tot{$gnames{$key}})."\n"; # need gene lengths here - this requires ensembl ID and transcript ID
			}
		}
	}
	close($out);
}

sub ReadFile{	# read in and sum the variants for each gene in a given chromosome
	my $infile      = $_[0];
	my $out         = $_[1];
	my $chr         = $_[2];
	my $g_length    = $_[3];
	my $g_length_tot= $_[4];
	my $to_include  = $_[5];
	my $out_het     = $_[6];
	my $found       = $_[7];
	my $gnames      = $_[8];
	my %count       = ();	# the number of variant-individual pairs for this gene
	my %icount	= ();	# the number of individuals carrying qualifying variants for this gene
	my %mis_mis     = ();   # number of missense only individuals for this gene
	my %mis_lof     = ();   # number of missense and lof individuals for this gene
        my %lof_lof     = ();   # number of lof_only individuals for this gene

	open(my $in,"results/chr_${chr}_variants_${infile}.counts") or die $!;	# open the results filr for this chromosome
	# read in gene names from the first line and save them to H
	my $line=<$in>;
	my @H=split(' ',$line);
	# now loop through each of the genes and initialise a hash element for icount for each of them, so that the hash element exists even if there isn't anyone with these
	for(my $i=2;$i<@H;$i++){
		my @F=split("_",$H[$i]);
		$icount{$H[$i]}=0;
		$mis_mis{$F[0]."_".$F[1]."_".$F[2]}=0;
                $mis_lof{$F[0]."_".$F[1]."_".$F[2]}=0;
                $lof_lof{$F[0]."_".$F[1]."_".$F[2]}=0;
	}
	# then loop through the file
	while(<$in>){
		chomp;
		my @F=split(' ');
        	my %person_mis  = ();
		my %person_lof  = ();
		if(exists $to_include->{$F[0]}){	# check if this person is in the list of individuals to include
			for(my $i=2;$i<@F;$i++){
				$count{$H[$i]}+=$F[$i];	# add variant-individual count to the hash
				if($F[$i]>0){	# and increment the individual count if the person has any number of variants
					$icount{$H[$i]}++;
					my @G=split("_",$H[$i]);
                                        if($G[3] eq "lof"){
                                                $person_lof{$G[0]."_".$G[1]."_".$G[2]}=$F[$i];
                                        }elsif($G[3] eq "missense"){
                                                $person_mis{$G[0]."_".$G[1]."_".$G[2]}=$F[$i];
					}
				}
			}
		}
		foreach my $key (keys %mis_mis){     # now count the number of lof/lof, lof/missense and missense/missense for this gene
                        if(exists $person_lof{$key} && exists $person_mis{$key}){
                                $mis_lof{$key}++;
                        }elsif(exists $person_lof{$key}){
                                $lof_lof{$key}++;
                        }elsif(exists $person_mis{$key}){
                                $mis_mis{$key}++;
			}
                }
	}
	close($in);
	foreach my $key (keys %count){
		my @F=split("_",$key);
		my @G=split("\\(",$F[0]);
		my @H=split("\\)",$G[1]);
		if(exists $g_length->{$G[0]."_".$F[1]}){
			print $out join("\t",$H[0],$G[0],$F[1],$F[2],$F[3],$count{$key},$icount{$key},$g_length->{$G[0]."_".$F[1]},$g_length_tot->{$G[0]."_".$F[1]})."\n";
		}else{
			print "WARNING: ".$key." not found in Ensembl gene list\n";
		}
		my @I=split(",",$H[0]);
		for(my $i=0;$i<@I;$i++){
			$$found{$I[$i]."_".$F[2]}=1;
		}
	}
	foreach my $key (keys %mis_mis){
                my @F=split("_",$key);
                my @G=split("\\(",$F[0]);
                my @H=split("\\)",$G[1]);
                if(exists $g_length->{$G[0]."_".$F[1]}){
                        print $out_het join("\t",$H[0],$G[0],$F[1],$F[2],$mis_mis{$key},$mis_lof{$key},$lof_lof{$key},$g_length->{$G[0]."_".$F[1]},$g_length_tot->{$G[0]."_".$F[1]})."\n";
                }else{
                        print "WARNING: ".$key." not found in Ensembl gene list\n";
                }
        }
}
