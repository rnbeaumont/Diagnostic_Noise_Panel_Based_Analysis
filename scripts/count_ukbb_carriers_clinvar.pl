#!/usr/bin/perl -w
# script to read in VEP reports, and list the number of G2P variants are carried by each UKBB individual
use strict;
use Getopt::Long;
use POSIX;
use Fcntl qw(SEEK_SET SEEK_CUR SEEK_END);

my $vep = undef;	# location of VEP reports
my $clinvar = undef;    # location of VEP report containing clinvar consequences
my $plink = undef;	# location of UKBB plink files
my $out = undef;	# output prefix
my $af_file_ukb = undef;	# file containing UKBB allele frequency report
GetOptions(
  "vep=s"     => \$vep,
  "plink=s"   => \$plink,
  "clinvar=s" => \$clinvar,
  "out=s"     => \$out,
  "freq_ukb=s" => \$af_file_ukb,
);

if ($vep && $plink && $out && $af_file_ukb && $clinvar) {
    &MainSub($vep, $plink, $out, $af_file_ukb, $clinvar);
}
else {
    print "\nconvert_vcf_to_snptest.pl\n";
    print "split up vcf containing imputation probabilities into nsnp long chunks\n";
    print "Usage:\n";
    print "        --vep       path to VEP report file\n";
    print "        --clinvar   file containing clinvar consequences from VEP\n";
    print "        --plink     prefix of plink files\n";
    print "        --freq_ukb  plink allele frequency report for UKBB\n";
    print "        --out       <output_file_prefix>\n\n";
}

sub MainSub{
    my $vep         = $_[0];
    my $plink       = $_[1];
    my $out         = $_[2];
    my $af_file_ukb = $_[3];
    my $clinvar     = $_[4];
    my $nsnpcount   = 0;
    my $nfilecount  = 1;
    my $header      = undef;
    my $file        = undef;

    # First read in the VEP reports and save the list of SNPs in a hash
    # Clinvar first
    print "Getting ClinVar consequences\n";
    my $clinvar_annotation=&ReadClinVarConsequence($clinvar);
    print "Reading VEP report\n";
    my ($snp_hash,$genes,$genes_hash,$snp_annotation)=&ReadVEPReport($vep);
    # read plink allele frequency
    print "Reading UKBB allele frequencies\n";
    my ($freq_ukb,$allele_ukb,$other_allele_ukb)=&ReadAF_ukb($af_file_ukb);
    # Next read in the Plink files, looking for the SNPs and counting the number which are carried by each person
    print "Reading plink files\n";
    &ReadPlink($plink,$snp_hash,$out,$genes,$genes_hash,$freq_ukb,$allele_ukb,$other_allele_ukb,$clinvar_annotation);
}

sub ReadPlink{
    my $plink         = $_[0];	# prefix for plink files to read
    my $snp_hash      = $_[1];	# hash containing the SNPs to look for
    my $out           = $_[2];	# output filename
    my $genes         = $_[3];	# array containing the names of the genes present
    my $genes_hash    = $_[4];	# hash containing the names of genes and whether they're biallelic or monoallelic
    my $freq_ukb      = $_[5];	# hash containing the allele frequencies in UKBB
    my $allele_ukb    = $_[6];	# hash containing the allele corresponding to the freq
    my $other_allele_ukb = $_[7];	# hash containing the allele corresponding to the freq
    my $snp_annotation = $_[8];	# hash of annotations in clinvar
    my $mono_freq     = 0.0001;	# allele frequency threshold for monoallelic genes
    my $biallelic_freq = 0.005;	# allele frequency threshold for biallelic genes
    my $bed           = undef;	# bed filehandle
    my $bim           = undef;	# bim filehandle
    my $fam           = undef;	# fam filehandle
    my $out_file      = undef;	# output filehandle
    my %gene_variant  = ();	# hash of variants in each gene for each person
    my %person_variant= ();	# hash of variant gene pairs for each person
    my $buffer        = undef;	# variable to hold the buffer from the bed file
    my $bits          = undef;	# the unpacked genotype data
    my @iid           = ();	# array of the IIDs from the fam file
    my $nfam          = undef;	# number of people in the fam file
    my $nbyte         = undef;	# number if bytes for each line
    my %found         = ();	# keep track of variants we've found

    my @C=("conflicting","benign","pathogenic","path-and-benign","uncertain");

    open($bed,$plink.".bed") or die $!;
    open($bim,$plink.".bim") or die $!;
    open($fam,$plink.".fam") or die $!;
    # read in the fam file
    while(<$fam>){
        chomp;
        my @F=split(' ');
        push @iid, $F[0];
    }
    $nfam=@iid;
    $nbyte=ceil($nfam/4);	# there are 4 individuals per byte (2 bits per individual) so we know how many bytes to skip per line
    binmode($bed);
    read($bed,$buffer,3);
    $bits=unpack("B*",$buffer);
    if($bits ne "011011000001101100000001"){
        die "ERROR: Plink .bed file in invalid format\n";
    }
    while(<$bim>){	# loop over the SNPs in the bim file
        chomp;
        my @F=split(' ');
	my @G=();	# array to store dosages in
	if(($. % 10000) == 0){	# print progress as we go
            print "SNP ".$.."\n";
	}
	# check if the SNP is in gnomad. If so, save the frequency, if not then set it to 0
	my $frq=0;
        if(exists $snp_hash->{$F[1]} && ($frq<$biallelic_freq || $frq>(1-$biallelic_freq))){	# if it's in the hash and the allele frequency is below the threshold then we need to process it
            $found{$F[1]}=1;
	    # check that the allele frequency corresponds to the alt allele here, otherwise flip it
	    my $flip=0;	# work out if it's coded to major allele and if so we'll need to flip it
	    my $gnomad=0;
	    $frq=$freq_ukb->{$F[1]};
	    if($gnomad==0){	# else, check which is the major allele in UKBB
		if($freq_ukb->{$F[1]}<0.5){
		    if($allele_ukb->{$F[1]} eq $F[4]){
			$flip=1;
		    }
		}else{
		    if($allele_ukb->{$F[1]} eq $F[5]){
			$flip=2;
		    }
		}
	    }
            for(my $i=0;$i<$nfam;$i+=4){	# read in the genotype data 4 individuals at a time
                read($bed,$buffer,1);
                $bits=unpack("B*",$buffer);	# unpack the SNP data to readable format
                my @H=split(//,$bits);	# and split it into an array
		for(my $j=7;$j>0;$j-=2){	# loop through the individuals in this byte, reading backward
                    if(!($H[$j]==1 && $H[$j-1]==0)){	# check if this individual is coded as missing
			my $dosage=undef;
			if($flip==1){
                            $dosage=(2-$H[$j]-$H[$j-1]);
		        }else{
                            $dosage=$H[$j]+$H[$j-1];
			    if($flip==2){
				    $frq=1-$frq;
			    }
			}
                        push @G, $dosage;
                    }else{
                        push @G, "0";
                    }
                }
            }
            # now we've got a vector of the genotypes of each person, so we can see if any carry the right dosage to be saved
	    my @I=split("\t",$snp_hash->{$F[1]});	# split out the gene names associated with this gene
	    my %variant_gene=();
	    for(my $i=0;$i<@I;$i++){	# put all the genes into a single hash (so we only count it once)
		$variant_gene{$I[$i]}=1;
	    }
	    # decide if the SNP is lof or missense
	    my $annot="none";
	    if(exists $snp_annotation->{$F[1]}){
		    $annot=$snp_annotation->{$F[1]};
	    }
	    # increment the annotation for this SNP appropriately
	    foreach my $key (keys %variant_gene){	# loop over the genes for this SNP
		my @H=split(',',$genes_hash->{$key});	# split out the required inheritance modes to test
		for(my $i=0;$i<@H;$i++){
                    if(($H[$i] eq "biallelic" && $frq<$biallelic_freq) || ($H[$i] eq "monoallelic" && $frq<$mono_freq)){	# check the allele frequency is OK, and if so count the carriers
                        for(my $j=0;$j<@iid;$j++){	# loop over the individuals and append the variants info for carriers to  an array, as well as counting the number if variants they have
                            $gene_variant{$key."_".$H[$i]."_".$iid[$j]."_".$annot}+=$G[$j];	# increment the number of variants in this gene for this person. Name is gene_interitance-mode_ID
                            $gene_variant{$key."_".$H[$i]."_".$iid[$j]}+=$G[$j];	# increment the number of variants in this gene for this person. Name is gene_interitance-mode_ID
	                    if($G[$j]>0){
                                $person_variant{$iid[$j]}.=$F[1]."_".$G[$j].",";	# save a list of all variants for this person
                                $person_variant{$key."_".$H[$i]."_".$iid[$j]."_".$annot}.=$F[1]."_".$G[$j].",";	# save a list of variants for this person for each gene
			    }
			}
                    }
		}
            }
        }else{
            seek($bed,$nbyte,SEEK_CUR);	# read in this SNP from the bed file to skip it
        }
    }
    # close the plink files
    close($fam);
    close($bim);
    close($bed);
    # and output the results
    print "Writing output files\n";

    ## counts of variants and "gene units" found for this person
    # output ID, chr total variants, and number of variants in each unit
    open($out_file,">",$out.".counts");
    # print the header
    print $out_file "ID\tchr_total_variants";
    for(my $i=0;$i<@$genes;$i++){
	my @F=split(",",$$genes_hash{$$genes[$i]});
	for(my $j=0;$j<@F;$j++){
	    for(my $c=0;$c<@C;$c++){
                print $out_file "\t".$$genes[$i]."_".$F[$j]."_".$C[$c];
	    }
	}
    }
    print $out_file "\n";
    # now loop through each person
    for(my $i=0;$i<@iid;$i++){
	print $out_file $iid[$i];	# output ID
        # first we need to count the total number of variants for this person
	my $total=0;
	if(exists  $person_variant{$iid[$i]}){	# if the person has any variants add the ones in valid units to a hash so we only count each variant once
            my %variants_valid=();
	    for(my $j=0;$j<@$genes;$j++){
		for(my $c=0;$c<@C;$c++){
		    if(exists $person_variant{$$genes[$j]."_monoallelic_".$iid[$i]."_".$C[$c]}){	# monoallelic count all variants
		        my @F=split(",",$person_variant{$$genes[$j]."_monoallelic_".$iid[$i]."_".$C[$c]});
		        for(my $k=0;$k<@F;$k++){
			    my @G=split("_",$F[$k]);
		            $variants_valid{$G[0]}=$G[1];
		        }
		    }
		}
		for(my $c=0;$c<@C;$c++){
		    if(exists $person_variant{$$genes[$j]."_biallelic_".$iid[$i]."_".$C[$c]}){	# biallelic so only count variants if there is more than 1
		        if($gene_variant{$$genes[$j]."_biallelic_".$iid[$i]}>1){
		            my @F=split(",",$person_variant{$$genes[$j]."_biallelic_".$iid[$i]."_".$C[$c]});
		            for(my $k=0;$k<@F;$k++){
			        my @G=split("_",$F[$k]);
		                $variants_valid{$G[0]}=$G[1];
			    }
		        }
		    }
		}
	    }
	    # now loop through the unique variants and sum their allele counts
	    foreach my $key (keys %variants_valid){
		$total+=$variants_valid{$key};
	    }
        }
        print $out_file "\t".$total;
	# now loop through each gene unit and output it
        for(my $j=0;$j<@$genes;$j++){
		my @F=split(",",$$genes_hash{$$genes[$j]});
		for(my $k=0;$k<@F;$k++){
			if($F[$k] eq "monoallelic"){
			    for(my $c=0;$c<@C;$c++){
            			if(exists $gene_variant{$$genes[$j]."_monoallelic_".$iid[$i]."_".$C[$c]}){
			        	print $out_file "\t".$gene_variant{$$genes[$j]."_monoallelic_".$iid[$i]."_".$C[$c]};
				}else{
					print $out_file "\t0";
				}
			    }
			}
			if($F[$k] eq "biallelic"){
			    for(my $c=0;$c<@C;$c++){
            			if(exists $gene_variant{$$genes[$j]."_biallelic_".$iid[$i]."_".$C[$c]}){
					if($gene_variant{$$genes[$j]."_biallelic_".$iid[$i]}>1){
		    				print $out_file "\t".$gene_variant{$$genes[$j]."_biallelic_".$iid[$i]."_".$C[$c]};
					}else{
						print $out_file "\t0";
					}
				}else{
		    			print $out_file "\t0";
	        		}
			    }
            		}
			if($F[$k] eq "imprinted"){
			    for(my $c=0;$c<@C;$c++){
            			if(exists $gene_variant{$$genes[$j]."_imprinted_".$iid[$i]."_".$C[$c]}){
					if($gene_variant{$$genes[$j]."_imprinted_".$iid[$i]}>1){
		    				print $out_file "\t".$gene_variant{$$genes[$j]."_imprinted_".$iid[$i]."_".$C[$c]};
					}else{
						print $out_file "\t0";
					}
				}else{
		    			print $out_file "\t0";
	        		}
	        	    }
            		}
			if($F[$k] eq "mosaic"){
			    for(my $c=0;$c<@C;$c++){
            			if(exists $gene_variant{$$genes[$j]."_mosaic_".$iid[$i]."_".$C[$c]}){
					if($gene_variant{$$genes[$j]."_mosaic_".$iid[$i]}>1){
		    				print $out_file "\t".$gene_variant{$$genes[$j]."_mosaic_".$iid[$i]."_".$C[$c]};
					}else{
						print $out_file "\t0";
					}
				}else{
		    			print $out_file "\t0";
	        		}
			    }
            		}
			if($F[$k] eq "digenic"){
			    for(my $c=0;$c<@C;$c++){
            			if(exists $gene_variant{$$genes[$j]."_digenic_".$iid[$i]."_".$C[$c]} && $F[$k]){
					if($gene_variant{$$genes[$j]."_digenic_".$iid[$i]}>1){
		    				print $out_file "\t".$gene_variant{$$genes[$j]."_digenic_".$iid[$i]."_".$C[$c]};
					}else{
		    				print $out_file "\t0";
					}
				}else{
		    			print $out_file "\t0";
	        		}
	        	    }
            		}
			if($F[$k] eq "uncertain"){
			    for(my $c=0;$c<@C;$c++){
            			if(exists $gene_variant{$$genes[$j]."_uncertain_".$iid[$i]."_".$C[$c]}){
					if($gene_variant{$$genes[$j]."_uncertain_".$iid[$i]}>1){
		    				print $out_file "\t".$gene_variant{$$genes[$j]."_uncertain_".$iid[$i]."_".$C[$c]};
					}else{
						print $out_file "\t0";
					}
				}else{
		    			print $out_file "\t0";
	        		}
	        	    }
            		}
	    	}
        }
	print $out_file "\n";
    }
    close($out_file);
    open($out_file,">",$out.".variant_info");
    print $out_file "ID\tvariants\n";
    foreach my $key (keys %gene_variant){
	if(exists $person_variant{$key}){
        	print $out_file join("\t",$key,$person_variant{$key})."\n";
	}
    }
    close($out_file);
    open($out_file,">",$out.".missing_variants");
    foreach my $key (keys %$snp_hash){
        if(!exists $found{$key}){
            print $out_file $key."\n";
        }
    }
    close($out_file);
}

sub ReadVEPReport{	# Read in the VEP report and save the results to a hash to return
    my $vep          = $_[0];
    my $file_line    = undef;
    my @data         = ();
    my @chr          = ();
    my %snp_hash     = ();
    my %snp_annotation = ();
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
                if($G[4] eq "-"){	# is's a deletion
                    my $pos=$G[1]-1;
                    my @I=split('',$G[3]);	# work out how many bases are deleted
                    my $ndel=@I;
                    if(!exists $snp_hash{$G[0].":".$pos.":D:".$ndel}){
                        $snp_hash{$G[0].":".$pos.":D:".$ndel}=$F[1]."_".$F[2];
                        $snp_annotation{$G[0].":".$pos.":D:".$ndel}=$G[6];
                    }else{
                        $snp_hash{$G[0].":".$pos.":D:".$ndel}.="\t".$F[1]."_".$F[2];
                        $snp_annotation{$G[0].":".$pos.":D:".$ndel}.=",".$G[6];
                    }
                }elsif($G[3] eq "-"){
                    my $pos=$G[1]-1;
                    my @I=split('',$G[3]);	# work out how many bases are inserted	
                    my $nins=@I;
	            if(!exists $snp_hash{$G[0].":".$pos.":I:".$nins}){
                        $snp_hash{$G[0].":".$pos.":I:".$nins}=$F[1]."_".$F[2];
                        $snp_annotation{$G[0].":".$pos.":I:".$nins}=$G[6];
			$nins--;
                        $snp_hash{$G[0].":".$G[1].":I:".$nins}=$F[1]."_".$F[2];
                        $snp_annotation{$G[0].":".$G[1].":I:".$nins}=$G[6];
                    }else{
                        $snp_hash{$G[0].":".$pos.":I:".$nins}.="\t".$F[1]."_".$F[2];
                        $snp_annotation{$G[0].":".$pos.":I:".$nins}.="\,".$G[6];
			$nins--;
                        $snp_hash{$G[0].":".$G[1].":I:".$nins}.="\t".$F[1]."_".$F[2];
                        $snp_annotation{$G[0].":".$G[1].":I:".$nins}.=",".$G[6];
                    }
                }else{	# else it's a SNP, so we can just save it
                    if(!exists $snp_hash{$G[0].":".$G[1].":".$G[3].":".$G[4]}){
                        $snp_hash{$G[0].":".$G[1].":".$G[3].":".$G[4]}=$F[1]."_".$F[2];	# so save whether it's het or hom to the hash
                        $snp_annotation{$G[0].":".$G[1].":".$G[3].":".$G[4]}=$G[6];
                    }else{
                        $snp_hash{$G[0].":".$G[1].":".$G[3].":".$G[4]}.="\t".$F[1]."_".$F[2];
                        $snp_annotation{$G[0].":".$G[1].":".$G[3].":".$G[4]}.=",".$G[6];
                    }
                }
            }
        }
    }
    close $snp_file;
    foreach my $key (keys %genes_hash){
        push @genes, $key;
    }
    return(\%snp_hash,\@genes,\%genes_hash,\%snp_annotation);
}

sub ReadAF{
    my $af_file = $_[0];	# plink allele frequency report
    my %freq    = ();	# hash of allele frequencies
    my %allele  = ();	# hash containing allele corresponding to frequency
    my %other_allele  = ();	# hash containing allele corresponding to frequency
    open(my $in,$af_file) or die $!;
    my $line=<$in>;
    while(<$in>){
        chomp;
	my @F=split(' ');
	# now check if it's an insertion or deletion
	my @H=split('',$F[2]);
	my @I=split('',$F[3]);
	if(@H>1){	# it's a deletion
		my $num=@H-1;
		$freq{$F[0].":".$F[1].":D:".$num}=$F[4];
		$allele{$F[0].":".$F[1].":D:".$num}=$F[3];
		$other_allele{$F[0].":".$F[1].":D:".$num}=$F[2];
	}elsif(@I>1){	# it's an insertion
		my $num=@I-1;
		$freq{$F[0].":".$F[1].":I:".$num}=$F[4];
		$allele{$F[0].":".$F[1].":I:".$num}=$F[3];
		$other_allele{$F[0].":".$F[1].":I:".$num}=$F[2];
	}else{	# is't a SNP
		$freq{$F[0].":".$F[1].":".$F[2].":".$F[3]}=$F[4];
		$allele{$F[0].":".$F[1].":".$F[2].":".$F[3]}=$F[3];
		$other_allele{$F[0].":".$F[1].":".$F[2].":".$F[3]}=$F[2];
		$freq{$F[0].":".$F[1].":".$F[3].":".$F[2]}=$F[4];	# save the name with allele order both ways around as we can't be sure which way UKBB have ordered them in the name
		$allele{$F[0].":".$F[1].":".$F[3].":".$F[2]}=$F[3];
		$other_allele{$F[0].":".$F[1].":".$F[3].":".$F[2]}=$F[2];
	}
    }
    close($in);
    return(\%freq,\%allele,\%other_allele);
}

sub ReadAF_ukb{
    my $af_file = $_[0];	# plink allele frequency report
    my %freq    = ();	# hash of allele frequencies
    my %allele  = ();	# hash containing allele corresponding to frequency
    my %other_allele  = ();	# hash containing allele corresponding to frequency
    open(my $in,$af_file) or die $!;
    <$in>;
    while(<$in>){
        chomp;
	my @F=split(' ');
	$freq{$F[1]}=$F[4];
	$allele{$F[1]}=$F[3];
	$other_allele{$F[1]}=$F[2];
    }
    close($in);
    return(\%freq,\%allele,\%other_allele);
}

sub ReadClinVarConsequence{
	my $clinvar = $_[0];	# clinvar file (VEP report with ClinVar consequences in
	open(my $in,"zcat $clinvar | ") or die $!;
	my %clinvar_cons=();
	my $line=undef;
	do{
		$line=<$in>;
	}while($line !~ /#Uploaded_variation/);
	while(<$in>){
		chomp;
		my @F=split(' ');
		# split out annotations and see if CLIN_SIG exists
		my @G=split(";",$F[13]);
		my $clinsig=undef;
		for(my $i=0;$i<@G;$i++){
			if($G[$i] =~ /CLIN_SIG/){
				my @H=split("=",$G[$i]);
				my $cpath=0;
				my $cben=0;
				my $cconf=0;
				my $unc=0;
				my @I=split(",",$H[1]);
				for(my $j=0;$j<@I;$j++){
				    if($I[$j] eq "pathogenic" || $I[$j] eq "pathogenic/likely_pathogenic" || $I[$j] eq "likely_pathogenic"){
					$cpath=1;
				    }
				    if($I[$j] eq "conflicting_interpretations_of_pathogenicity"){
					    $cconf=1;
				    }
				    if($I[$j] eq "uncertain_significance"){
					    $unc=1;
				    }
				    if($I[$j] eq "likely_benign" || $I[$j] eq "benign/likely_benign" || $I[$j] eq "benign"){
					    $cben=1;
				    }
				}
				if($unc){
					$clinsig="uncertain";
				}elsif($cconf){
					$clinsig="conflicting";
				}elsif($cben && $cpath){
				    $clinsig="path-and-benign"
				}elsif($cben==1){
				    $clinsig="benign";
			        }else{
				    $clinsig="pathogenic";
				}
			}
		}
		if($clinsig){
			if(!exists $clinvar_cons{$F[0]}){
				$clinvar_cons{$F[0]}=$clinsig;
			}else{
				if($clinvar_cons{$F[0]} eq "pathogenic"){
					$clinvar_cons{$F[0]}=$clinsig;
				}elsif($clinvar_cons{$F[0]} eq "benign" && $clinsig eq "conflicting"){
					$clinvar_cons{$F[0]}=$clinsig;
				}
			}
		}
	}
	close($in);
	return(\%clinvar_cons);
}
