#!/usr/bin/perl -w
# run the script: time perl 00.script/00.select.ncRNA.pl 01.data/00.PriorData 01.data/00.PriorData/total.RNA.fasta Arabidopsis

use strict;

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $tgtfile = shift @ARGV;
my $sp = shift @ARGV;

## run the script
opendir(SRC,"$srcfolder") || die"Cannot open file:$!";
my @sams = sort(grep(/total.RNA.part[0-9]+.fasta/, readdir(SRC)));
closedir SRC;
open(TGT, ">>$tgtfile");

my $count=0;
foreach my $sam (@sams){
	if($sam eq "total.RNA.part1.fasta"){
		system("cat $srcfolder/total.RNA.part1.fasta > $tgtfile");
	}else{
		open(SAM, "$srcfolder/$sam");
		my $mark = 0;
		
		while(my $line = <SAM>){
			chomp $line;
			if($line =~ /^>/){
				if($line =~ /$sp/){
					$count++;
					print TGT ">$count\n";
					$mark = 1;
				}else{
					$mark = 0;
				}
			}else{
				if($mark == 1){
					print TGT "$line\n";
				}
			}
		}
		
		close SAM;
	}
}

close TGT; 