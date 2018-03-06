#!/usr/bin/perl -w
# run the script: time perl 00.script/04.folder.ncRNA.pl 01.data/01.Fastq/SRA174408 01.data/03.CutAdapter 01.data/05.Blat 01.data/06.ncRNA Sapelo 10

use strict;
system("echo 'Running 04.folder.ncRNA.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $reffolder = shift @ARGV;
my $cutfolder = shift @ARGV;
my $srcfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $mode = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;

my @temp = split(/\//, $reffolder);
my $sample = pop @temp;
my $prev = "00.script/03.blat.script/$sample";
my $dir = "00.script/04.ncRNA.script1/$sample";
=cut
## check if previous step has succesfully finished
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^\w+/, readdir(CHK)));
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("blat.$sample.$chk.sh.o*");
		my @stdout = glob("blat.$sample.$chk.sh.e*");
		my @log = glob("$prev/$chk.done.*");
		if(!($stderr[0] and $stdout[0] and $log[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' blat.$sample.$chk.sh.e* >> $prev/summary.error.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "$prev/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$srcfolder/$chk/$chk.bls") and (!(-s "$srcfolder/$chk/$chk.R1.bls") or !(-s "$srcfolder/$chk/$chk.R2.bls"))){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					#system("rm -f blat.$chk.*");
					#system("rm -f $prev/$chk.done.log");
					system("echo 'Resubmitting the job: blat.$chk.sh' >> job.monitor.txt");
					#system("qsub 00.script/03.blat.script/blat.$chk.sh");
						#last; # There are some jobs failed
				}
				else{
					$i++;
				}
			}
			if($i == scalar @chks){
				system("echo 'All jobs have been finished successfully' >> job.monitor.txt"); # error file is empty
				last;
			}
		}else{
			die "ERROR: something went wrong in previous steps\n";	
		}
	}
	sleep $sleeptime;
}
closedir CHK;
=cut
## start running the script
opendir(SRC,"$srcfolder") || die"Cannot open file:$!";
my @sams = sort(grep(/\w+/, readdir(SRC)));
closedir SRC;

system("mv blat.$sample.* $prev/");
system("rm -rf $dir");
system("mkdir -p $dir");

foreach my $sam(@sams)
{
	my $shell = "$dir/ncRNA.$sample.$sam.sh";
	open(SHL,">$shell") or die "Cannot write $shell: $!";
	
	if($platform eq "sapelo"){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N ncRNA.$sample.$sam.sh\n";
		print SHL "#PBS -l nodes=1:ppn=1:AMD\n";
		print SHL "#PBS -l walltime=12:00:00\n";
		print SHL "#PBS -l mem=30gb\n\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	
	print SHL "mkdir -p $tgtfolder/$sam\n";
	
	if($mode eq "single-end"){
		#print SHL "time perl 00.script/04.clrna.pl $srcfolder/$sam/$sam.bls $cutfolder/$sam/$sam.fastq $tgtfolder/$sam/$sam.fastq\n";
        print SHL "echo \$((`wc -l < $tgtfolder/$sam/$sam.fastq` / 4 )) > $tgtfolder/$sam/$sam.count.txt\n";
	}elsif($mode eq "paired-end"){
		#print SHL "time perl 00.script/04.clrna.pl $srcfolder/$sam/$sam.R1.bls $cutfolder/$sam/$sam.R1.fastq $tgtfolder/$sam/$sam.R1.fastq\n";
		#print SHL "time perl 00.script/04.clrna.pl $srcfolder/$sam/$sam.R2.bls $cutfolder/$sam/$sam.R2.fastq $tgtfolder/$sam/$sam.R2.fastq\n";
        print SHL "echo \$((`wc -l < $tgtfolder/$sam/$sam.R1.fastq` / 4 )) > $tgtfolder/$sam/$sam.R1.count.txt\n";
        print SHL "echo \$((`wc -l < $tgtfolder/$sam/$sam.R2.fastq` / 4 )) > $tgtfolder/$sam/$sam.R2.count.txt\n";
	}
	
	print SHL "touch $dir/$sam.done.log\n";
    close SHL;
	system("chmod u+x $shell");
	if($platform eq "sapelo"){
    	system("qsub $shell");
	}elsif($platform eq "zcluster"){
		system("qsub -q rcc-30d $shell");
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
 
}

system("echo 'Finished 04.folder.ncRNA.pl!' >> job.monitor.txt");
