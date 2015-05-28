#!/usr/bin/perl                                         
#Columbia test program - Ariella Sasson

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper qw(Dumper);

$0 =~ s-.*/--g;

my $title="";
my ($fna_file,$qual_file,$output_file);

GetOptions( 'f|fna=s'  => \$fna_file,
	    'q|qual=s'  => \$qual_file,
	    'o|output=s' => \$output_file
	   );

print "\n";

my @assign;
my $acnt;
@assign=&AssignInputs();
$acnt=@assign;

my $primer_seq="CGCCGTTTCCCAGTAGGTCTC";
my $adaptor_seq="ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG"; 

if ($acnt > 0){
    my $trash = join('',@assign);
    print "\n$trash";
    usage();
}

print "Inputs \n";
print "fna File: $fna_file\n"; 
print "Quality File: $qual_file\n";
print "Output File prefix: $output_file \n";
print "Primer Sequence: $primer_seq\n";
print "Primer Sequence length: ".length($primer_seq)."\n";
print "Adaptor Sequence: $adaptor_seq\n";
print "Adapter Sequence length: ".length($adaptor_seq)."\n";

print "\n"."Starting Analysis ....\n";

{
#load in files
my @readArray = &read_multi_fna($fna_file);
my @qualArray = &read_multi_qual($qual_file);
my $read_cnt = @readArray;

#generate the files with the primer and adapter seq
open (PSEQ, ">blast_seq.txt");
print PSEQ ">primer\n$primer_seq\n";
print PSEQ ">adapter\n$adaptor_seq\n";
close(PSEQ);

#blast commands - generate db and blast primer and adapter seq
`formatdb -i $fna_file -p F`;
`blastall -p blastn -d test.fna -i blast_seq.txt -m 8 -o blast_seq_results.txt`;

#read in blast m8 output and create a hash
print "\nBlast output file in m8 format: blast_seq_results.txt\n";
my %blast_hoa = &read_blast_m8("blast_seq_results.txt");

#output files
open(FNA100,">".$output_file."_rawRead_gt100_adapterPrimer_Trimmed.fna");
open(QUAL100,">".$output_file."_rawRead_gt100_adapterPrimer_Trimmed.qual");

open(FNA20,">".$output_file."_rawRead_qavg_20_adapterPrimer_Trimmed.fna");
open(QUAL20,">".$output_file."_rawRead_qavg_20_adapterPrimer_Trimmed.qual");

open(BLAST,">".$output_file."_blast_adapterPrimer_Location.txt");
print "FNA file for raw reads w length gt 100: ".$output_file."_rawRead_gt100_adapterPrimer_Trimmed.fna\n";
print "Qual file for raw reads w length gt 100: ".$output_file."_rawRead_gt100_adapterPrimer_Trimmed.qual\n";
print "FNA file for raw reads w quality average gt 20: ".$output_file."_rawRead_qavg_20_adapterPrimer_Trimmed.fna\n";
print "Qual file for raw reads w quality average gt 20: ".$output_file."_rawRead_qavg_20_adapterPrimer_Trimmed.qual\n";
print "Tab delimited primer and adapter locations: ".$output_file."_blast_adapterPrimer_Location.txt\n";

#counters
my $j=0;
my $raw_avg20_cnt=0;
my $raw_length_cnt=0;
my $read_w_primer=0;
my $read_w_adapter=0;
my $reads_w_both=0;

#analysis
while($j<$read_cnt){
    my $seq=$readArray[$j][1];
    my @qual_split=split(" ",$qualArray[$j][1]);
    my $seq_len=length($seq);
    my $q_len=@qual_split;
    my $sum = 0;
    $sum += $_ for (@qual_split);
    my $avg=$sum/$q_len;
    my $cntr2=0;
    my $cntr1=0;
    my $adj_seq=$seq;
    my @adj_qual_split=@qual_split;

#checks raw quality >=20
    if($avg>=20){
	$cntr2++;
	$raw_avg20_cnt++;
    }
    
#checks raw length >=100
    if($seq_len>=100){
	$cntr1++;
	$raw_length_cnt++;
    }

    my @holda=split(" ",$readArray[$j][0]);
    my $name = shift(@holda);
    my $ind_p=0;
    my $ind_a=0;

#checks hash if there was a blast hit
    if(exists($blast_hoa{$name})){
	my @phold = split("\n", $blast_hoa{$name});
	my @blast_pos;
	my $z=0;

	foreach(@phold){
	    my @blast_h=split("\t",$_);
	    my $per_id=$blast_h[2];
	    if ($per_id < 90){
		next;
	    }

	    my @pos_hold=($blast_h[8],$blast_h[9]);
	    my @sorted=sort @pos_hold;
	    push(@blast_pos,[@sorted]);

	    $z++;
	    print BLAST "$blast_h[1]\t$sorted[0]\t$sorted[1]\n";

	    if ($blast_h[0] eq "primer"){
		$read_w_primer++;
		$ind_p=1;
	    }
	    if ($blast_h[0] eq "adapter"){
                $read_w_adapter++;
                $ind_a=1;
            }
	}#foreach(@phold){

	if ($ind_a==1 & $ind_p==1){
	    $reads_w_both++;
	}

	my @sorted_pos = sort { $b->[1] <=> $a->[1]} @blast_pos;

#trims read
	foreach my $row (@sorted_pos) {
	    my $start_pos = @$row[0];
	    my $end_pos = @$row[1];
	    my $length_pa=$end_pos-$start_pos+1;
	    substr($adj_seq, $start_pos-1, $length_pa)="";
	    splice(@adj_qual_split, $start_pos-1, $length_pa);

	}
#	print Dumper \@sorted_pos;

    }#if(exists($blast_hoa{$name})){

    if($cntr1 > 0){
	print FNA100 ">$name\n$adj_seq\n";
	print QUAL100 ">$name\n@adj_qual_split\n";
      }

    if($cntr2 > 0){
        print FNA20 ">$name\n$adj_seq\n";
        print QUAL20 ">$name\n@adj_qual_split\n";
    }

    $j++;
}


print "\n\nNumber of reads in the datase: $j\n";
print "Number of raw reads with avg quality > 20: $raw_avg20_cnt\n";
print "Number of raw reads with length > 100: $raw_length_cnt\n";
print "Number of reads with primer: $read_w_primer\n";
print "Number of reads with adapter: $read_w_adapter\n";
print "Number of reads with primer and adapter: $reads_w_both\n";

exit;
}

###################################
########SUBROUTINES################
###################################

sub AssignInputs{
    my @err;

    if ((!defined($fna_file)) || ($fna_file eq '')) {
        push(@err, "Error 1: FNA file is not defined\n");
    }

    if ((!defined($qual_file)) || ($qual_file eq '')) {
	push(@err, "Error 2: Qual File is not defined\n");
    }
	    
    if ((!defined($output_file)) || ($output_file eq '')) {
	push(@err, "Error 3: Output File is not defined\n");
    }
    return @err;
}

sub read_multi_fna{
    my $inFile = shift; 
    my @reads_array;
    my $i=0;
    open (IN, "$inFile");

    local $/ = ">";  #slurp the fasta file

    my $garbage = <IN>; # Discard the ">" at the begining of the file

    while ( my $record = <IN> ) {
	chomp $record; 
	my ($name, @seqLines) = split /\n/, $record;
	
	my $sequence = join('',@seqLines); # Concatenates all elements of the @seqLines array into a single string.

	$reads_array[$i][0]=$name;
	$reads_array[$i][1]=$sequence;
#	print "$i\n$name\n$sequence\n##################\n";
	$i++;

    }#end while
    return @reads_array;
}

sub read_multi_qual{
    my $inFile = shift;
    my @qual_array;
    my $i=0;
    open (IN, "$inFile");

    local $/ = ">";  #slurp the fasta file                                      

    my $garbage = <IN>; # Discard the ">" at the begining of the file           

    while ( my $record = <IN> ) {
        chomp $record;
	my ($name, @qualLines) = split /\n/, $record;

        my $quality = join(' ',@qualLines); # Concatenates all elements of the @seqLines array into a single string.

        $qual_array[$i][0]=$name;
        $qual_array[$i][1]=$quality;
#        print "$i\n$name\n$quality\n##################\n";
        $i++;

    }#end while
    return @qual_array;
}

sub read_blast_m8{
    my $inFile = shift;
    my @blast_array;
    open (IN, "$inFile");
    my %masterHash;

    my $i=0;
    my $j=0;
    while ( my $record = <IN> ) {
        chomp $record;
	$i++;
	my @hold = split("\t",$record);

	
	if(exists($masterHash{$hold[1]})){
            my $hold1 = $masterHash{$hold[1]};
	    my @holda = split("\t",$hold1);
	    $hold1 = $hold1."\n".$record;

            $masterHash{$hold[1]}=$hold1;

	} else {
            $masterHash{$hold[1]}=$record;
	    $j++;
        }

    }#end while

    return %masterHash;
}

sub usage {

    print "\nusage: $0 \n[-f file.fna \n -q file.qual \n -o output_prefix ) \n] \n\n";
    print "This script takes an fna and quality file and blasts it against a known adapter and primer seq (embedded in the program) and the prints out the blast output, the fna and qual files with the adapter and primers trimmed where the raw reads have a length gt 100, the fna and qual files with the adapter and primers trimmed where the raw reads have a average quality gt 20, and the output of the read names with the start and stop locations in the read of the blast output.  This file only looks at blast matches of greater than 90% identity. \n\n";

    print "-f|fna - fna file - required \n\n";
    print "-q|qual - qual file for the fna file  - required \n\n";
    print "-o|output_prefix - output prefix name  - required \n\n";


    print "example: $0 -f file.fna -q file.qual -o out_o \n\n";

    exit;

}
