#!/usr/bin/perl
use warnings;
use strict;
my $dir="./";
my @files=glob($dir."*/*.gz");
my %variants=();
for(my $i=0;$i<@files;$i++){
    open(my $in, "zcat $files[$i] | ") or die $!;
    my $line=<$in>;
    while($line !~ /#CHROM/){
        $line=<$in>;
    }
    while(<$in>){
        chomp;
        my @F=split(' ');
        my @a=sort($F[3],$F[4]);
        $variants{$F[0].":".$F[1].":".$F[2].":".$a[0].":".$a[1]}=join("\t",$F[0],$F[1],$F[2],$F[3],$F[4],"42.5","PASS",".","GT:GQ:MIN_DP","0/0:30:30","0/1:30:30","1/1:30:30");
    }
    close($in);
}
open(my $out,">","DDD_variant_list_unordered");
print $out "##fileformat=VCFv4.2\n";
print $out join("\t","#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","ID_1","ID_2","ID_3\n");
foreach my $key (keys %variants){
    print $out $variants{$key}."\n";
}
close($out);
