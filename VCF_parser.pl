# Written April 4, 2012 by the Bioinformatics Core at the Children's Hospital of Philadelphia
#modified July 23 2014 for SNPeff 3.6
# Transforms any vcf into a tab delimited file
# WARNING:  this program assumes a well formed VCF.  It uses the header to get the column descriptions for the info field and the Sample descriptors.

#!/usr/bin/perl                                         
                                                                           
use strict;
use warnings;
use Getopt::Long;

$0 =~ s-.*/--g;

my $title="";
my ($in_vcf,$outputFile);

GetOptions( 'v|vcf=s'  => \$in_vcf,
	   );

print "\n";

my @assign;
my $acnt;
@assign=&AssignInputs();
$acnt=@assign;

$outputFile = $in_vcf;
$outputFile =~ s/.vcf$/_tab.txt/;


if ($acnt > 0){
    my $trash = join('',@assign);
    print "\n$trash";
    usage();
}

print "Inputs \n";
print "VCF File: $in_vcf\n"; 
print "Output SNP File Name: $outputFile \n";

print "\n"."Starting Analysis ....\n";

{
calcCoverageStats($in_vcf,$outputFile);

exit;
}

###################################
########SUBROUTINES################
###################################

sub AssignInputs{
    my @err;

    if ((!defined($in_vcf)) || ($in_vcf eq '')) {
        push(@err, "Error 1: VCF file is not defined\n");
    }

    return @err;
}

sub calcCoverageStats {
    my($in,$out) = @_;

    unless(open(TAB, "+> $out")){
        print "can't open $out to write to - please check your location.\n";
        return;
    }

    open(VCF, "$in") or die "can't open up your VCF file, $in . please check your directory location ";
 
   print "we have found the  file: $in \n";

    my %info_fields;
    my @info_headers;
    my @info_values;
    my @gt_headers;
    my @main_head;
    my @samples;
    my $sample_ind=0;

    push(@gt_headers,"Genotype Value (GT)");
    push(@gt_headers,"Genotype (GT)");
    push(@gt_headers,"Genotype Description (GT)");
    push(@gt_headers,"Genotype Quality (GQ)");
    push(@gt_headers,"Quality Filtered Depth (DP)");
    push(@gt_headers,"Raw Depth Ref (AD)");
    push(@gt_headers,"Raw Depth Alt (AD)");
    push(@gt_headers,"Phred Scaled Likelihood (PL)");
#    push(@gt_headers,"Phred Scaled Likelihood Het (PL)");
#    push(@gt_headers,"Phred Scaled Likelihood Hom Alt (PL)");

    my $eff_name="";
    my $eff_cnt=0;
    my $snpEff_check=0;

    my $gt_ind=1;

#parse Header Row
    while(<VCF>){ # it will either be one chromosome or all
        chomp;
	my $hold = $_;

	if(/^\#\#/ && !eof(VCF)){
	    if (/SnpEffVersion/){
		$snpEff_check=1;
	    }
	    if (/^\#\#INFO/){
		my @hold=split("<",$_);
		my @hold2=split(",",$hold[1]);
		foreach (@hold2){
		    if (/ID/){
			my @holdID = split("=",$_);
			if ($holdID[1] =~ m/EFF/ && $hold !~ m/Not provided/){
			    $snpEff_check=2;
			    my @form=split("'",$hold);
			    my @effect=split('\(',$form[1]);
			    my @eff2=split('\)',$effect[1]);
			    my @eff3=split('\[',$eff2[0]);
			    my @effDet=split('\|',$eff3[0]);
			    $effect[0]=~s/\s+//g;
			    $eff_name="SNPEff.".$effect[0];
			    $eff_cnt++;
			    foreach (@effDet){
				s/\s+//g;
				$eff_cnt++;
				$eff_name=$eff_name."\t"."SNPEff.".$_;
			    }
			    $eff_name=$eff_name."\t"."SNPEff.[".$eff3[1];
			    $eff_cnt++;
			}
		        $info_fields{$holdID[1]} = "" unless exists $info_fields{$holdID[1]};
		    }
		}
	    }
	    next;
        }elsif(/^\#CHROM/ && !eof(VCF)) {
	    my @hold = split("\t",$_);
	    my $pos=0;
	    my $len = @hold;
	    if ($len == 8){
		$gt_ind=0;
	    }
	    for (my $i=0;$i<$len;$i++){
		my $col = $hold[$i];
		if ( $col !~ m/INFO/){
		    if ( $col !~ m/FORMAT/){
			push(@main_head,$col);
		    }
		}
		if($col=~/FORMAT/){
		    $pos=$i+1;
		    last;
		}
	    }
	    if ($gt_ind==1){
		
		for (my $i=$pos;$i<$len;$i++){
		    my $sample= $hold[$i];
		    push(@samples,$sample);
		}
		my $slen = @samples;
		if ($slen >= 1){
		    $sample_ind=1;
		}
		last;
	    }else{
		last;
	    }
	}
    }
    print "\nHeader Parsing Complete....\n";

#if there are no sample names gives it one - will assume only a single sample in this case
    if ($sample_ind == 0 && $gt_ind ==1 ){
	push(@samples,"sample");
    }

    my $effind=0;
    my $effval="";
    my $eff_description_title = "\'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] )\'";

    if ($snpEff_check==1){
	$effind=1;
	$effval="EFF";;
	my @form=split("'",$eff_description_title);
	my @effect=split('\(',$form[1]);
	my @eff2=split('\)',$effect[1]);
	my @eff3=split('\[',$eff2[0]);
	my @effDet=split('\|',$eff3[0]);
	$effect[0]=~s/\s+//g;
	$eff_name="SNPEff.".$effect[0];
	$eff_cnt++;
	foreach (@effDet){
	    s/\s+//g;
	    $eff_cnt++;
	    $eff_name=$eff_name."\t"."SNPEff.".$_;
	}
	$eff_name=$eff_name."\t"."SNPEff.[".$eff3[1];
#	print " $eff_name \n @info_headers\n";
	$eff_cnt++;
    }

#    print "$eff_description_title \n";

    foreach my $key (sort keys %info_fields) {

	if ($key eq "EFF"){
	    $effind=1;
	    $effval=$key;
	    next;
	}
	push(@info_headers,$key);
	push(@info_values,"");
    }

    if ($effind==1){
	push(@info_headers,$effval);
    }
    
    
    my @GT_title;
    if ($gt_ind ==1){
	foreach my $s (@samples){
	    foreach my $g (@gt_headers){
		push(@GT_title,$s.".".$g);
	    }
	}
    }
    
    my $a;
    if ($effind==1){
	$a=pop (@info_headers);
	push (@info_headers,$eff_name);
    }

    my $hI = join("\t",@info_headers);
    my $hG = join("\t",@GT_title);
    my $hM = join("\t",@main_head);

    if ($effind==1){
	pop @info_headers;
	push @info_headers,$a;
    }
    
	print TAB $hM."\t".$hI."\t".$hG."\n";


    print "\nStarting to parse VCF entries....\n\n";
    my $vcf_line_cnt=0;
    my $tab_line_cnt=0;
    
    while(<VCF>){
#	print "$_\n";
	$vcf_line_cnt++;
	my ($mv,$gt,@info_values)=&readVCFLine($_,$gt_ind,\@info_headers);
	my $iv = "";

	if ($effind==1 && !(m/EFF\=/)){
	    $tab_line_cnt++;
	    $iv=join("\t",@info_values);
	    for (my $n=1;$n< $eff_cnt;$n++){
		$iv = $iv."\t N\/A";
	    }
	    print TAB $mv."\t".$iv.$gt."\n";
	}elsif ($effind==1){
	    my $eff_p_val = pop(@info_values);
	    $iv=join("\t",@info_values);
	    my @val_eff_p = &parse_Eff($eff_p_val,$eff_name);
	    foreach my $p(@val_eff_p){
		$tab_line_cnt++;
		print TAB $mv."\t".$iv."\t".$p.$gt."\n";
	    }
	}else{
	    $tab_line_cnt++;
	    $iv=join("\t",@info_values);
	    print TAB $mv."\t".$iv.$gt."\n";
	}

    } # while (<FH>)
    print "\nTotal Number of vcf entries parsed: $vcf_line_cnt\n";
    print "Total Number of tabbed entries written: $tab_line_cnt\n\n";

    close(TAB);
    close(VCF);
}

sub parse_Eff{
    my ($eff_val,$eff_hdr) = @_;
    my @hold_name = split("\t",$eff_hdr);
    my $cnt = @hold_name;
    my @hold_eff_val=split(",",$eff_val);
    my $cnt2=@hold_eff_val;
    my $value="";
    my @name_arr;

    foreach my $entry (@hold_eff_val){
	my $ecnt=0;
	my @effect=split('\(',$entry);
	my @eff2=split('\)',$effect[1]);
	my @eff3=split('\[',$eff2[0]);
	my $c = @eff3;
	my @effDet=split('\|',$eff3[0]);
	my $d= @effDet;

	$value=$effect[0];
	$ecnt++;
	foreach (@effDet){
	    $ecnt++;
	    $value=$value."\t".$_;
	}

	if ($ecnt < $cnt-1){
	    for(my $i = $ecnt+1;$i<$cnt;$i++){
		$ecnt++;
		$value=$value."\t"."";
	    }
	}
	if ($c==1){
	    $value=$value."\t"."[\|\|\]";
	    $ecnt++;
	}else{
	    $ecnt++;
	    $value=$value."\t"."[".$eff3[1];
	}
	push( @name_arr,$value);
    }
    return(@name_arr);
}

sub getGT{
    my ($genoLabels,$rb,$sb,$geno) = @_;
    my $typeAlt="N/A";
    my $typeAltDesc="N/A";
    my $description="N/A";
    my $rawDP = "N/A";
    my $rawAD_ref="N/A";
    my $rawAD_alt="N/A";
    my $GQ="N/A";
    my $PL_val="N/A";
#    my $PL_homA="N/A";
#    my $PL_het="N/A";

    if ($geno !~ /\.\/\./){
	if ($genoLabels ne ""){
	    my @name = split(':',$genoLabels);
	    my @type = split(':',$geno);
	    
	    for my $i (0 .. $#name){
		if ($name[$i] =~ /GT/){
		    $typeAlt = $type[$i];
		    my @gttype=split("/",$typeAlt);
		    if ($gttype[0]==$gttype[1]){
			$description = "Hom";
		    }else{
			$description = "Het";
		    }
		    my @sb_v=split(",",$sb);
		    unshift(@sb_v,$rb);

		    $typeAltDesc = "$sb_v[$gttype[0]]"."/"."$sb_v[$gttype[1]]"; 
#		    print "$geno $typeAlt -- $sb $rb -- @sb_v -- $typeAltDesc\n";

		}elsif($name[$i] =~ /AD/){
		    my @cov_ad = split(",",$type[$i]);
		    my $len = @cov_ad;
		    
		    $rawAD_ref =$cov_ad[0];
#		    if ($len >2){
			shift @cov_ad;
			$rawAD_alt = join(",",@cov_ad);
#		    }else{
#			$rawAD_alt =$cov_ad[1];
#		    }
		}elsif($name[$i] =~ /DP/){
		    $rawDP=$type[$i];
		}elsif($name[$i] =~ /GQ/){
		    $GQ=$type[$i];
		}elsif($name[$i] =~ /PL/){

		    $PL_val=$type[$i];
#		    my @phred_PL = split(",",$type[$i]);
                    
                    #Actual VCF values as seen in the file.  I fyou want the calculated value use the commented out lines below
#                    $PL_homR=$phred_PL[0];
#                    $PL_homA=$phred_PL[2];
#		    $PL_het=$phred_PL[1];

##		    $PL_homR=sprintf("%.2f", 10**((-1)*($phred_PL[0]/10)));
##		    $PL_homA=sprintf("%.2f", 10**((-1)*($phred_PL[2]/10)));
##		    $PL_het=sprintf("%.2f", 10**((-1)*($phred_PL[1]/10)));
		}
	    }
	}
    }else{
	$typeAlt="./.";
	$typeAltDesc = "."."/".".";
	$description = "No Call Possible";
    }
#    my $holdgt = $typeAlt."\t".$typeAltDesc."\t".$description."\t".$GQ."\t".$rawDP."\t".$rawAD_ref."\t".$rawAD_alt."\t".$PL_homR."\t".$PL_het."\t".$PL_homA;
    my $holdgt = $typeAlt."\t".$typeAltDesc."\t".$description."\t".$GQ."\t".$rawDP."\t".$rawAD_ref."\t".$rawAD_alt."\t".$PL_val;
    return $holdgt;

}
sub readVCFLine {
    my ($pileupLine,$gtInd,$Iheader,$Ival) = @_;

    if (!$pileupLine) {
        return ("","","","");
    }
    chomp($pileupLine);

    my @row = split(/\t/, $pileupLine); #split the line

    my $len=@row;
    my $len_head=scalar(@{$Iheader});

    my @info = split(";",$row[7]);
    my @info_val=("N/A") x $len_head;

    my $genoLabels=$row[8];
    my $GT_sample_val="";

my ($currPileupChrom,$currPileupPos,$id,$refBase,$seqBase,$qual,$filter) = ($row[0],$row[1],$row[2],$row[3],$row[4],$row[5],$row[6]);

    if ($gtInd ==1){
	for(my $i=9;$i<$len;$i++){
	    my $gt_val = &getGT($genoLabels,$refBase,$seqBase,$row[$i]);
	    $GT_sample_val = $GT_sample_val."\t".$gt_val;
	}
    }
#    my ($currPileupChrom,$currPileupPos,$id,$refBase,$seqBase,$qual,$filter) = ($row[0],$row[1],$row[2],$row[3],$row[4],$row[5],$row[6]);

    my $main_val = $currPileupChrom."\t".$currPileupPos."\t".$id."\t".$refBase."\t".$seqBase."\t".$qual."\t".$filter;

    for (my $m=0;$m<$len_head;$m++){
	my $item = @$Iheader[$m];
	foreach my $vcf (@info){
	    if ($vcf =~ /^$item=/){
		my @holdpv = split("=",$vcf);
		$info_val[$m]=$holdpv[1];
	    }
	}
    }

    return ($main_val,$GT_sample_val,@info_val);

}





sub usage {

    print "\nusage: $0 \n[-v raw vcf file]\n\n";
    print "This script takes the variant vcf file and and transforms it to a tab delimited file.  If the annotation is done with SnpEff, the vcf EFF entry is also parsed and tab delimited.  Merged VCFs is supported.\nOutput file changes the .vcf ending to _tab.txt\n";

    print "-v|vcf - vcf file - required \n\n";

    print "example: $0 -v A_var_raw.vcf \n\n";

    exit;

}
