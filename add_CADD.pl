#!/usr/bin/perl
use strict;
use warnings;

$0 =~ s-.*/--g;

# sorry, using globals
my $tabix="/nas/is1/bin/tabix";

use DBI;
if ($#ARGV != 2 ) {
    print "assumes tabix CADD  files\n";
    print "usage: perl add_CADD.pl input.vcf output.vcf cadd.txt.gz\n";
    exit;
}

my $file=$ARGV[0];
my $outfile = $ARGV[1];
my $cadd=$ARGV[2];


open (VCFF, $file);
open (OVCFF, ">",$outfile);
#open (CADD,$cadd);

my $totalcount=0;
my $cadd_counter=0;

while (my $vcfline=<VCFF>) {
    chomp $vcfline;
    
    #skipping the comment lines
    if($vcfline =~ m/#/){
	if($vcfline =~ m/^##reference=/i){
	    print OVCFF "$vcfline\n";
	    print OVCFF "##INFO=<ID=CADD_RAW,Number=1,Type=Float,Description=\"CADD Combined Annotation Dependent Depletion - Raw Score - from $cadd\">\n";
	    print OVCFF "##INFO=<ID=CADD_SCALED,Number=1,Type=Float,Description=\"CADD Combined Annotation Dependent Depletion - Scaled Score - from $cadd\">\n";
	    next;
	}
	print OVCFF "$vcfline\n";
	next;	
    }
    
    my @line = split("\t",$vcfline);
    my $mychr = $line[0];
    my $chr = $mychr;
    $chr=~s/chr//g;
    $totalcount++;
    
    my $pos = $line[1];
    my $ref = $line[3];
    my $alt = $line[4];
    my $linefull = $vcfline;

    my $ref_len = length $ref;
    
    if ($ref_len > 1){
	print OVCFF "$vcfline\n";
	next;
    }

    my @alt_array = split(",",$alt);
    my $len_alt_array = @alt_array;
    
#    print "-----------------\n $line[0] $line[1] $line[2] $line[3] $line[4] $line[5] $line[6] \n $line[7]\n $alt --- $len_alt_array $alt @alt_array\n";

    my $location = $chr.":".$pos."-".$pos;
    my $cadd_raw="";
    my $cadd_scaled="";

    my $cadd_results = `$tabix $cadd $location`;


 #   print "+++++++++++++++++++++++\nafter:\n Sift: $sift_results \n PP2 $pp2_results \ncadd\n$cadd_results\n";

    my @cadd_results_array=split("\n",$cadd_results);
    
    foreach my $alt_res (@alt_array){
	my $alt_var_len=length $alt_res;
	if($alt_var_len > 1){
	    next;
	}
	
	foreach my $result (@cadd_results_array) {
	    my @split_result=split("\t",$result);
	    my $r_chrom=$split_result[0];
	    my $r_pos=$split_result[1];
	    my $r_ref=$split_result[2];
	    my $r_alt=$split_result[3];
	    my $r_raw=$split_result[4];
	    my $r_scaled=$split_result[5];
	    
	    if ($r_alt eq $alt_res){
		$cadd_raw=$cadd_raw.$r_raw;
		$cadd_scaled=$cadd_scaled.$r_scaled;
		$cadd_counter++;
		last;
	    }
	    
	}#foreach my $result (@cadd_results_array) {
	
	$cadd_raw=$cadd_raw.",";
	$cadd_scaled=$cadd_scaled.",";
	
    } #foreach my $alt_res (@alt_array){
    
#    print "Before $chr $pos $ref $alt \t$cadd_raw\t$cadd_scaled\t$sift_score\t$pp2_score\n";

    if (($cadd_raw ne "")&($cadd_raw ne ",")){
	$cadd_raw = substr($cadd_raw, 0, -1);
        $cadd_scaled= substr($cadd_scaled, 0, -1);
	$linefull=~ s/\tGT:/;CADD_RAW=$cadd_raw;CADD_SCALED=$cadd_scaled\tGT:/;
    }

    print OVCFF "$linefull\n";

} #while (my $vcfline=<VCFF>)

print "Total variants:$totalcount\n";
print "CADD Annotated variants:$cadd_counter\n";
#print "SIFT Annotated variants:$sift_counter\n";
#print "PP2 Annotated variants:$pp2_counter\n";

close (VCFF);
close (OVCFF);
exit;
