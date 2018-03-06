#!/bin/perl -w

#use strict;

if (@ARGV<3) {
	die "$0:Not enough arguments\n";
}

my $frblat=$ARGV[0];
#my $fwmatch=$ARGV[1];
my $frseq=$ARGV[1];
my $fwseq=$ARGV[2];
my $line;
my $seq;
my @temp;
my @rna1;
my @rna2;
my $id;
my $name1;
my $name2;

my $count1=0;
my $count2=0;
my $mark=0;
my $i=0;
my $test=0;
my $size = 50000000;

#test the global arguments.
my $len=@ARGV;
print $len,"\n";
#print $fwseq,"\n";
print "Blat file:$ARGV[0];\nFastq file=$ARGV[1];\nCleaned file=$ARGV[2]\n";

#The first round, to pick up the pvalue >=90.

open(FRB,"$frblat")||die "Can not open the blat file: $!.\n";
#open(FWM,">$fwmatch")||die "Can not create the highly matched file: $!.\n";
#open(FRM,"$fwmatch")||die "Can not open the highly matched file: $!.\n";

my $k = 0;
my $t = 0;
my $ncrna_num = 0;
my $fq_num=0;
my $fq_out_num=0;

while(1){
	#Build a hash to keep the sequence ID
	my %hash=();
	$t ++;
	while ($line = <FRB>) {
		$k ++;
		chomp $line;
		@temp=split(/\t/,$line);
		$name1=$temp[0];
	#	$name2=$temp[1];	
		my @seq_id=split('_',$name1);	
		my $id=$seq_id[0];
		my $len=$seq_id[1];
		my $len_hit=$temp[2]*$temp[3]/100;

		if ($len_hit/$len>=0.9) {
			$hash{$id}=1;
			$ncrna_num ++;
		}#first if
		
		if($k % $size == 0){
			last;
		}
	}#first foreach

	#print FWM "The total number of matches are: $count1.\n";
	print "$k\n";
	print "Processing part $t\n";
	#close(FWM);

	#The second round, to remove the matched sequence.
	$count1=0;

	if($t == 1){
		open(FRS,"$frseq")||die "Can not open the RNAseq file: $!.\n";
		open(TGT, ">$fwseq.$t");
	}else{
		my $prevt = $t - 1;
		if($t > 2){
			my $rm = $t - 2;
			system("rm $fwseq.$rm");
		}
		open(FRS,"$fwseq.$prevt")||die "Can not open the RNAseq file: $!.\n";
		open(TGT, ">$fwseq.$t");
	}

	my $fq_id='';
	my $fq_out='';
	my $fq_len=0;
	$fq_num=0;
	$fq_out_num=0;

	while  ($seq=<FRS>) {
			if ($seq=~ m/^\@/) {
				$fq_num++;
				@temp=split(/\s+/,$seq);
				$fq_id=substr ($temp[0],1);
				$fq_out=$seq;
				$seq=<FRS>;
				$fq_len=length($seq);
				$fq_out.=$seq;
				for (1..2){
					$seq = <FRS>; 
					$fq_out.=$seq;
				}#for
			}#if
			
			## check if ID is in the list
			if((not exists $hash{$fq_id}) && ($fq_len>=25)){
				$fq_out_num++;
				print TGT $fq_out;
			}
	}
	
	close TGT;
	
	if(!$line){
		last;
	}
}

system("mv $fwseq.$t $fwseq");
system("rm $fwseq.*");
close(FRB);

print "total reads: $fq_num\n ncRNA reads: $ncrna_num\n left reads: $fq_out_num\n";
        
        
			 
			 