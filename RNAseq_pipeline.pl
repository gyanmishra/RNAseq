#!/usr/bin/perl

#### Gyan Prakash Mishra, j12mishra@gmail.com 
#### RNA-seq data analysis using Tophat(v2.1.0) and Cufflink(v2.2.1)
#### User can provide no. of thread for running all the command ,rest uses default parameters
#### Usage: $ perl RNAseq_pipeline.pl <sample file > <Path File>
#### Example: perl RNAseq_pipeline.pl sample.txt file_path.txt


use strict ;
use warnings;
use Data::Dumper;
use POSIX;
use Cwd qw();

my $sampleinfo     =$ARGV[0];
my $file_path      =$ARGV[1];

open(FH, "$sampleinfo") or die;
open(FH1, "$file_path") or die; 

my @sample = <FH> ;
my @file_path= <FH1> ;
my $fastq1="";
my $fastq2="";
my $bowtieindex ="";
my $reference_gtf ="";
my $sample ="";
my $stimulation ="";
my $noofrep ="";
my $replicate ="";
my $output	="";
my $arg ="";
my @bam_file ;
my $thread ="";
#my $replicate1 ="";
my $pwd = getcwd();

#print "$pwd\n";


#################################################################################################################################
foreach my $line(@file_path) {
	chomp($line);
	unless ($line =~ /^#/) {
		my @temp = split("=", $line);
		$temp[0] =~ s/\r|\n|\s//g;
		$temp[1] =~ s/\r|\n|\s//g;
		chomp($temp[0],$temp[1]);
		#print "$temp[0]\t$temp[1]\n";

		if($temp[0] eq "bowtie_index")
		{
			$bowtieindex = $temp[1];
			
		}
		if($temp[0] eq "reference_gtf")
		{
			$reference_gtf = $temp[1];
		}

	}
}
#print "$bowtieindex\t$reference_gtf\n";

################################################################################################################################

################################################################################################################################
#
#							Run tophat2					
#							###########
print " Does your sample have paired end reads : (Y/N):";
$arg =<STDIN>;
print " How many replicates you have for each sample : ";
$noofrep = <STDIN>;
print "How many threads you want to use to run tophat:";
$thread =<STDIN>;
foreach my $line1(@sample)
{
	chomp($line1);
	my @csv = split /\t/, $line1;
	$stimulation = $csv[0];
        $sample = $csv[1];
	
	if($arg eq "Y" )
	{
		$fastq1  = $csv[2];
		$fastq2	= $csv[3];
		$replicate =$csv[4];
			
		print "$csv[4]\n";
		$output = $sample."_".$replicate;
		#tophat($fastq1,$fastq2,$thread,$bowtieindex, $reference_gtf);
		my $run_tophat_pe = "tophat -p $thread-o tophat_out/$stimulation/$output -G $reference_gtf $bowtieindex $fastq1 $fastq2";
	        push @bam_file ,"tophat_out/$stimulation/$output";
		print "Running TOPHAT\n";
		print "$run_tophat_pe\n";
		system($run_tophat_pe);
		
		#print "\n$fastq1\n$fastq2\n";
	}
	else 
	{
		$fastq1  = $csv[2];
		$replicate =$csv[3];
		$output = $sample."_".$replicate;
		my $run_tophat_se = "tophat -p $thread -o tophat_out/$stimulation/$output -G $reference_gtf $bowtieindex $fastq1";
                push @bam_file ,"tophat_out/$stimulation/$output";
		print "Running TOPHAT\n";
		print "$run_tophat_se\n";
		system($run_tophat_se);
	}
}

################################################################################################################################

################################################################################################################################
#
#							Run Cufflink
#							############

while(<FH>)
{
	chomp($_);
	
	my $cufflink ="cufflinks -o $stimulation/$sample -G $reference_gtf tophat_out/$stimulation/$sample/accepted_hits.bam";
 	print "$cufflink\n";	
	system($cufflink)
}

#################################################################################################################################

#################################################################################################################################
#
#							Run Cuffmerge
#							#############
open(IN, "$sampleinfo") or die;
open(OUT, ">","assembly_list.txt") or die"Couldn't open: $!" ;

while(<IN>)
{
	chomp($_);
        my @csv1 = split /\t/, $_;
	$stimulation = $csv1[0];
	my $rep = $csv1[3];
	my $output1 =$sample."_".$rep ;

	print OUT "tophat_out/$stimulation/$output1/transcripts.gtf\n";
}

my $merged = "Control"."_"."Treatment"."_"."merged" ;
my $cuffmerge ="cuffmerge -o $merged/$sample -g $reference_gtf assembly_list.txt";
print "$cuffmerge\n";
system($cuffmerge)

#################################################################################################################################
#################################################################################################################################
##
##               					Run Cuffdiff
#							############
if($noofrep == 2)
{
	my $cuffmerge ="cuffdiff -p $thread -o $merged $merged/$sample/merged.gtf -L control,$sample 
	$bam_file[0]/accepted_hits.bam,$bam_file[1]/accepted_hits.bam 
	$bam_file[2]/accepted_hits.bam,$bam_file[3]/accepted_hits.bam";
	system($cuffmerge);
}
if($noofrep == 3)
{
        my $cuffmerge ="cuffdiff -p $thread -o $merged $merged/$sample/merged.gtf -L control,$sample 
	$bam_file[0]/accepted_hits.bam,$bam_file[1]/accepted_hits.bam,$bam_file[2]/accepted_hits.bam 
	$bam_file[3]/accepted_hits.bam,$bam_file[4]/accepted_hits.bam,$bam_file[5]/accepted_hits.bam";
        system($cuffmerge);
}

##################################################################################################################################
#
#		Filter raw gene_exp.diff from cuffdiff output for list of Significantly Differentiated genes
#		############################################################################################			#

#open(IN1, "$merged/gene_ex.diff") or die "Coudn't open this file, There could be some problem in previous steps";
print "Cuffdiff Analysis is complete !! \n";
print "Extract Differentially expressed genes";
print "Please proveide q value[e.g 0.01,0.05 etc] and Log2(Fold change)[e.g 1,0.5 etc] cutoff to get Differentially expressed gene\n";
print "q value:";
my $qvalue=<STDIN>;
print "log2(Fold change)";
my $lfc=<STDIN>;

if($qvalue == 0.05 && $lfc == 1))
{
	`grep -P "OK|gene_id" $merged/$sample/gene_exp.diff | sort -k 13n,13n | perl -ne '@data=split("\t", $_); if ($data[12]<=0.05){print;}' >DE_genes.txt`;
	`awk -F"\t" ' {if($10 >1) print }' DE_genes.txt >2_fold_Upregulated_gene`;
	`awk -F"\t" ' {if($10 < -1) print }' DE_genes.txt >2_fold_Downregulated_gene`;

}
if($qvalue == 0.01 && $lfc == 1))
{
        `grep -P "OK|gene_id" $merged/$sample/gene_exp.diff | sort -k 13n,13n | perl -ne '@data=split("\t", $_); if ($data[12]<=0.01){print;}' >DE_genes.txt`;
        `awk -F"\t" ' {if($10 >1) print  }' DE_genes.txt >2_fold_Upregulated_gene`;
        `awk -F"\t" ' {if($10 < -1) print }' DE_genes.txt >2_fold_Downregulated_gene`;
}
if($qvalue == 0.05 && $lfc == 0.5))
{
        `grep -P "OK|gene_id" $merged/$sample/gene_exp.diff | sort -k 13n,13n | perl -ne '@data=split("\t", $_); if ($data[12]<=0.05){print;}' >DE_genes.txt`;
        `awk -F"\t" ' {if($10 >0.5) print }' DE_genes.txt >1_fold_Upregulated_gene`;
        `awk -F"\t" ' {if($10 < -0.5) print }' DE_genes.txt >1_fold_Downregulated_gene`;

}
if($qvalue == 0.01 && $lfc == 0.5))
{
        `grep -P "OK|gene_id" $merged/$sample/gene_exp.diff | sort -k 13n,13n | perl -ne '@data=split("\t", $_); if ($data[12]<=0.01){print;}' >DE_genes.txt`;
        `awk -F"\t" ' {if($10 >0.5) print  }' DE_genes.txt >1_fold_Upregulated_gene`;
        `awk -F"\t" ' {if($10 < -0.5) print }' DE_genes.txt >1_fold_Downregulated_gene`;
}




	
	

