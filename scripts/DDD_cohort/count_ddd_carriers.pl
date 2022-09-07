#!/usr/bin/perl -w
# script to read in VEP reports, and list the number of G2P variants are carried by each UKBB individual
use strict;
use Getopt::Long;
use POSIX;
use Fcntl qw(SEEK_SET SEEK_CUR SEEK_END);

my $vep = undef;	# location of VEP reports
my $plink = undef;	# location of UKBB plink files
my $out = undef;	# output prefix
my $af_file = undef;	# file containing gnomad allele frequency data
my $ddd_af = undef;	# optional file containing DDD allele frequenecies
GetOptions(
  "vep=s"     => \$vep,
  "plink=s"   => \$plink,
  "out=s"     => \$out,
  "freq=s"    => \$af_file,
  "ddd_freq=s"=> \$ddd_af,
);

if($vep && $plink && $out) {
    &MainSub($vep, $plink, $out);
}
else {
    print "\nconvert_vcf_to_snptest.pl\n";
    print "split up vcf containing imputation probabilities into nsnp long chunks\n";
    print "Usage:\n";
    print "        --vep       path to VEP report file\n";
    print "        --plink     prefix of plink files\n";
    print "        --out       <output_file_prefix>\n\n";
}

sub MainSub{
    my $vep         = $_[0];
    my $plink       = $_[1];
    my $out         = $_[2];
    my $nsnpcount   = 0;
    my $nfilecount  = 1;
    my $header      = undef;
    my $file        = undef;

    # First read in the VEP reports and save the list of SNPs in a hash
    print "Reading VEP report\n";
    my ($snp_hash,$genes,$genes_hash)=&ReadVEPReport($vep);
    # read plink allele frequency
    # Next read in the Plink files, looking for the SNPs and counting the number which are carried by each person
    print "Reading plink files\n";
    &ReadVCFFiles($plink,$snp_hash,$out,$genes,$genes_hash);
}

sub ReadVCFFiles{
    my $plink         = $_[0];	# prefix for plink files to read
    my $snp_hash      = $_[1];	# hash containing the SNPs to look for
    my $out           = $_[2];	# output filename
    my $genes         = $_[3];	# array containing the names of the genes present
    my $genes_hash    = $_[4];	# hash containing the names of genes and whether they're biallelic or monoallelic
    my $mono_freq     = 0.0001;	# allele frequency threshold for monoallelic genes
    my $biallelic_freq = 0.005;	# allele frequency threshold for biallelic genes
    my $bed           = undef;	# bed filehandle
    my $bim           = undef;	# bim filehandle
    my $fam           = undef;	# fam filehandle
    my $out_file      = undef;	# output filehandle
    my $buffer        = undef;	# variable to hold the buffer from the bed file
    my $bits          = undef;	# the unpacked genotype data
    my @iid           = ();	# array of the IIDs from the fam file
    my $nfam          = undef;	# number of people in the fam file
    my $nbyte         = undef;	# number if bytes for each line
    my %found         = ();	# keep track of variants we've found

    # setup output files
    ## counts of variants and "gene units" found for this person
    # output ID, chr total variants, and number of variants in each unit
    open(my $out_count,">",$out.".counts");
    # print the header
    print $out_count "ID\tchr_total_variants";
    for(my $i=0;$i<@$genes;$i++){
	my @F=split(",",$$genes_hash{$$genes[$i]});
	for(my $j=0;$j<@F;$j++){
            print $out_count "\t".$$genes[$i]."_".$F[$j];
	}
    }
    print $out_count "\n";
    open(my $out_var,">",$out.".variant_info");
    print $out_var "ID\tvariants\n";

    my @files=glob($plink."/*/*.gz");	# get list of DDD files
    for(my $i=0;$i<@files;$i++){	# loop through the individuals
        print $files[$i]."\n";
        my %gene_variant  = ();	# hash of variants in each gene for each person
        my %person_variant= ();	# hash of variant gene pairs for each person
        open(my $in,"zcat $files[$i] | ") or die $!;
        # skip header
        my $line=<$in>;
        while($line !~ /#CHROM/){
            $line=<$in>;
        }
        chomp $line;
        my @head=split(' ',$line);	# save the header with the ID of the individual
	# now process all the SNPs
        while(<$in>){
            chomp;
            my @F=split(' ');
            # build SNP ID
            # first check if it's an indel
	    my @A2=split(",",$F[4]);
	    for(my $a2=0;$a2<@A2;$a2++){
                my $num_1=split('',$F[3]);
                my $num_2=split('',$A2[$a2]);
                my $id=$F[0].":".$F[1].":".$F[3].":".$A2[$a2];
                if($num_1==1 && $num_2==1 && $A2[$a2] ne "*" && $A2[$a2] ne "-"){	# it's s sensibly named SNP
                    $id=$F[0].":".$F[1].":".$F[3].":".$A2[$a2];
                }elsif($A2[$a2] eq "<DEL>" || $A2[$a2] eq "-" || $A2[$a2] eq "*"){	# it's a deletion
                    $id=$F[0].":".$F[1].":D:".$num_1;
                }elsif($F[3] eq "<DEL>" || $F[3] eq "-"){	# insertion
                    $id=$F[0].":".$F[1].":I:".$num_2;
                }elsif($A2[$a2] eq "<DUP>"){	# duplication (insertion)
                    my $pos=$F[1]+$num_1;
                    $id=$F[0].":".$pos.":I:".$num_2;
                }elsif($num_1==1 && $num_2>$num_1){	# it's an insertion and the position needs incrementing
                    my $num=$num_2-$num_1;
                    my $pos=$F[1]+1;
                    $id=$F[0].":".$pos.":I:".$num;
                }elsif($num_1>$num_2 && $num_2==1){	# it's an deletion and the position needs incrementing
                    my $num=$num_1-$num_2;
                    my $pos=$F[1]+1;
                    $id=$F[0].":".$pos.":D:".$num;
                }elsif($num_1==$num_2){	# it's a SNP but multiple alleles so trim these and correct the position if needed
                    my @A=split('',$F[3]);
                    my @B=split('',$A2[$a2]);
                    my $pos=-1;
                    for(my $i=0;$i<@A;$i++){
                        if($A[$i] ne $B[$i]){
                            if($pos<0){
                                $pos=$i;
                            }else{
                                #print $F[0].":".$F[1].":".$F[3].":".$A2[$a2]."\n";
                            }
                        }
                        my $pnew=$F[1]+$pos;
                        $id=$F[0].":".$pnew.":".$A[$pos].":".$B[$pos];
                    }
                }elsif($num_1>$num_2){	# it's a deletion with multiple alleles, so trim off the remaining ones and update the position
                    my $len=$num_1-$num_2;
                    my $pos=$F[1]+$num_2;
                    $id=$F[0].":".$pos.":D:".$len;
                }else{	# else the second must be longer so it's an insertion
                    my $len=$num_2-$num_1;
                    my $pos=$F[1]+$num_1;
                    $id=$F[0].":".$pos.":I:".$len;
                }
                my $frq=0;	# initialise freq variable
                my $valid=1;	# does it pass QC? Assume yes until we find otherwise
                my $dddaf=undef;	# INFO field DDD_AF
                my @A=split(";",$F[7]);	# look for AF and quality metrics
                for(my $i=0;$i<@A;$i++){
                    my @B=split("=",$A[$i]);
                    if($B[0] eq "DDD_AF"){
		        my @Z=split(",",$B[1]);
		        $dddaf=$Z[$a2];
			if($dddaf eq "."){
			    $dddaf=0;
			}
                    }
                }
                if($dddaf){
                    $frq=$dddaf;
                }
                if((exists $snp_hash->{$id}) && ($frq<$biallelic_freq || $frq>(1-$biallelic_freq))){	# check if it's potentially one we want to include
                    $found{$id}=1;	# keep track of variants we've found
                    # now that we know we want to process this SNP, split out the gene names associated with this gene
                    my @I=();
                    if(exists $snp_hash->{$id}){
                        @I=split("\t",$snp_hash->{$id});
                    }
                    my %variant_gene=();	# because the gene names could appear twice, hash them to make sure we don't double count
                    for(my $i=0;$i<@I;$i++){
                        $variant_gene{$I[$i]}=1;
                    }
                    foreach my $key (keys %variant_gene){       # loop over the genes for this SNP
	                my @H=();
                        if(exists $genes_hash->{$key}){
                            @H=split(',',$genes_hash->{$key});	# split out the required inheritance modes to test
                        }
                        for(my $i=0;$i<@H;$i++){
                            if(($H[$i] eq "biallelic" && $frq<$biallelic_freq) || ($H[$i] eq "monoallelic" && $frq<$mono_freq)){	# check the allele frequency is OK, and if so count the carriers
                                my @Z=split(":",$F[9]);
                                my @G=split("/",$Z[0]);	# split out the genotypes of the child
                                my $dose=0;	# add together variant dosages
				if(@G==2){
				    if($G[0] == $a2+1){
					$dose++; # add together variant dosages
				    }
				    if($G[1] == $a2+1){
					$dose++; # add together variant dosages
				    }
				}else{
				    #$dose=2;
				}
                                $gene_variant{$key."_".$H[$i]}+=$dose;	# increment the number of variants in this gene for this person. Name is gene_interitance-mode
                                if($dose>0){
                                    $person_variant{$head[9]}.=$F[0].":".$F[1]."_".$dose.",";	# save a list of all variants for this person
                                    $person_variant{$key."_".$H[$i]}.=$F[0].":".$F[1]."_".$dose.",";	# save a list of variants for this person for each gene
		                }
                            }
		        }
		    }
                }
            }
        }
        close($in);
        # now print out this person
        print $out_count $head[9];	# output ID
        # first we need to count the total number of variants for this person
	my $total=0;
        my %variants_valid=();
	for(my $j=0;$j<@$genes;$j++){
	    if(exists $person_variant{$$genes[$j]."_monoallelic"}){	# monoallelic count all variants
	        my @F=split(",",$person_variant{$$genes[$j]."_monoallelic"});
	        for(my $k=0;$k<@F;$k++){
	            my @G=split("_",$F[$k]);
		    $variants_valid{$G[0]}=$G[1];
		}
            }
	    if(exists $person_variant{$$genes[$j]."_biallelic"}){	# biallelic so only count variants if there is more than 1
	        if($gene_variant{$$genes[$j]."_biallelic"}>1){
	            my @F=split(",",$person_variant{$$genes[$j]."_biallelic"});
	            for(my $k=0;$k<@F;$k++){
	                my @G=split("_",$F[$k]);
		        $variants_valid{$G[0]}=$G[1];
		    }
		}
            }
	}
        foreach my $key (keys %variants_valid){
            $total+=$variants_valid{$key};
        }
        print $out_count "\t".$total;
	# now loop through each gene unit and output it
        for(my $j=0;$j<@$genes;$j++){
            my @F=split(",",$$genes_hash{$$genes[$j]});
	    for(my $k=0;$k<@F;$k++){
		if($F[$k] eq "monoallelic"){
                    if(exists $gene_variant{$$genes[$j]."_monoallelic"}){
			print $out_count "\t".$gene_variant{$$genes[$j]."_monoallelic"};
		    }else{
			print $out_count "\t0";
		    }
		}
		if($F[$k] eq "biallelic"){
            	    if(exists $gene_variant{$$genes[$j]."_biallelic"}){
			if($gene_variant{$$genes[$j]."_biallelic"}>1){
		            print $out_count "\t".$gene_variant{$$genes[$j]."_biallelic"};
			}else{
			    print $out_count "\t0";
			}
		    }else{
		    	print $out_count "\t0";
	            }
            	}
		if($F[$k] eq "imprinted"){
            	    if(exists $gene_variant{$$genes[$j]."_imprinted"}){
			if($gene_variant{$$genes[$j]."_imprinted"}>1){
		            print $out_count "\t".$gene_variant{$$genes[$j]."_imprinted"};
			}else{
			    print $out_count "\t0";
			}
		    }else{
		    	print $out_count "\t0";
	            }
            	}
		if($F[$k] eq "mosaic"){
            	    if(exists $gene_variant{$$genes[$j]."_mosaic"}){
			if($gene_variant{$$genes[$j]."_mosaic"}>1){
		            print $out_count "\t".$gene_variant{$$genes[$j]."_mosaic"};
			}else{
			    print $out_count "\t0";
			}
		    }else{
		    	print $out_count "\t0";
	            }
            	}
		if($F[$k] eq "digenic"){
            	    if(exists $gene_variant{$$genes[$j]."_digenic"} && $F[$k]){
			if($gene_variant{$$genes[$j]."_digenic"}>1){
		            print $out_count "\t".$gene_variant{$$genes[$j]."_digenic"};
			}else{
		    	    print $out_count "\t0";
			}
		    }else{
		    	print $out_count "\t0";
	            }
            	}
		if($F[$k] eq "uncertain"){
            	    if(exists $gene_variant{$$genes[$j]."_uncertain"}){
			if($gene_variant{$$genes[$j]."_uncertain"}>1){
		    	    print $out_count "\t".$gene_variant{$$genes[$j]."_uncertain"};
			}else{
			    print $out_count "\t0";
			}
		    }else{
		    	print $out_count "\t0";
	            }
            	}
	    }
        }
        print $out_count "\n";
        print $out_var $head[9]."\t".$person_variant{$head[9]}."\n";
    }
    close($out_count);
    close($out_var);
}

sub ReadVEPReport{	# Read in the VEP report and save the results to a hash to return
    my $vep          = $_[0];
    my $file_line    = undef;
    my @data         = ();
    my @chr          = ();
    my %snp_hash     = ();
    my @genes        = ();
    my %genes_hash   = ();

    open(my $snp_file, $vep) or die "Unable to open file $vep\n";	# open the file
    while($file_line=<$snp_file>){	# loop through the lines
        chomp $file_line;
        my @F=split("\t",$file_line);	# split out the fields
        if($F[3] eq "is_canonical"){	# check whether it's labelled as in the canonical transcript
		#	    push @genes, $F[1]."_".$F[2];
            my @J=split('=',$F[5]);	# split REQ field
            my $biallelic = $J[1];
	    $genes_hash{$F[1]."_".$F[2]}=$biallelic;
	    my @H=split(";",$F[6]);	# split up the variants in this gene and loop through them
	    for(my $i=0;$i<@H;$i++){
                my @G=split(":",$H[$i]);	# split out the fields in the variant annotation field
                if($G[0] eq "X"){	# recode X and Y chr codes so they match VEP output
                    $G[0]="23";
                }elsif($G[0] eq "Y"){
                    $G[0]="24";
                }
                # now save the variant info
                # First we need to check if it's an indel. If so we'll need to fix the positions so they match the plink file
		my @A=split('',$G[3]);
		my @B=split('',$G[4]);
		my $num1=@A;
		my $num2=@B;
                if($G[4] eq "-"){	# is's a deletion
                    my $pos=$G[1];
                    my @I=split('',$G[3]);	# work out how many bases are deleted
                    my $ndel=@I;
                    if(!exists $snp_hash{$G[0].":".$pos.":D:".$ndel}){
                        $snp_hash{$G[0].":".$pos.":D:".$ndel}=$F[1]."_".$F[2];
                    }else{
                        $snp_hash{$G[0].":".$pos.":D:".$ndel}.="\t".$F[1]."_".$F[2];
                    }
                }elsif($G[3] eq "-"){
                    my $pos=$G[1];
                    my @I=split('',$G[3]);	# work out how many bases are inserted	
                    my $nins=@I;
	            if(!exists $snp_hash{$G[0].":".$pos.":I:".$nins}){
                        $snp_hash{$G[0].":".$pos.":I:".$nins}=$F[1]."_".$F[2];
			$nins--;
                        $snp_hash{$G[0].":".$G[1].":I:".$nins}=$F[1]."_".$F[2];
                    }else{
                        $snp_hash{$G[0].":".$pos.":I:".$nins}.="\t".$F[1]."_".$F[2];
			$nins--;
                        $snp_hash{$G[0].":".$G[1].":I:".$nins}.="\t".$F[1]."_".$F[2];
                    }
	        }elsif($num2>$num1){	# Insertion
		    my $nins=$num2-$num1;
		    my $pos=$G[1]+$num1;
		    if(!exists $snp_hash{$G[0].":".$pos.":I:".$nins}){
		        $snp_hash{$G[0].":".$pos.":I:".$nins}=$F[1]."_".$F[2];
		    }else{
		        $snp_hash{$G[0].":".$pos.":I:".$nins}.="\t".$F[1]."_".$F[2];
		    }
	        }elsif($num1>$num2){	# Deletion
		    my $ndel=$num1-$num2;
		    my $pos=$G[1]+$num2;
		    if(!exists $snp_hash{$G[0].":".$pos.":D:".$ndel}){
		        $snp_hash{$G[0].":".$pos.":D:".$ndel}=$F[1]."_".$F[2];
		    }else{
		        $snp_hash{$G[0].":".$pos.":D:".$ndel}.="\t".$F[1]."_".$F[2];
		    }
	        }elsif($num1==$num2 && $num1>1){
		    my @A=split('',$G[3]);
		    my @B=split('',$G[4]);
		    my $pos=-1;
		    for(my $i=0;$i<@A;$i++){
			if($A[$i] ne $B[$i]){
			    if($pos<0){
				$pos=$i;
			    }
			}
		    }
		    my $pnew=$G[1]+$pos;
		    $snp_hash{$G[0].":".$pnew.":".$A[$pos].":".$B[$pos]}=$F[1]."_".$F[2];
                }else{	# else it's a SNP, so we can just save it
                    if(!exists $snp_hash{$G[0].":".$G[1].":".$G[3].":".$G[4]}){
                        $snp_hash{$G[0].":".$G[1].":".$G[3].":".$G[4]}=$F[1]."_".$F[2];	# so save whether it's het or hom to the hash
                        $snp_hash{$G[0].":".$G[1].":".$G[4].":".$G[3]}=$F[1]."_".$F[2];	# so save whether it's het or hom to the hash
                    }else{
                        $snp_hash{$G[0].":".$G[1].":".$G[3].":".$G[4]}.="\t".$F[1]."_".$F[2];
                        $snp_hash{$G[0].":".$G[1].":".$G[4].":".$G[3]}.="\t".$F[1]."_".$F[2];
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
