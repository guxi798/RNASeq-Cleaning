#!/usr/bin/perl -w
# run the script: time perl 00.script/05.trimmomatic.pl 01.data/01.Fastq/SRA174408 01.data/06.ncRNA 01.data/07.cleanRNA Sapelo 10 01.data/00.PriorData/SRA174408

use strict;
system("echo 'Running 05.trimmomatic.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $reffolder = shift @ARGV;
my $srcfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $mode = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my $special = shift @ARGV;

my @temp = split(/\//, $reffolder);
my $sample = pop @temp;
my $prev = "00.script/04.ncRNA.script/$sample";
my $dir = "00.script/05.trimmomatic.script/$sample";

## check if previous step has succesfully finished
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^\w+/, readdir(CHK)));
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("ncRNA.$sample.$chk.sh.o*");
		my @stdout = glob("ncRNA.$sample.$chk.sh.e*");
		my @log = glob("$prev/$chk.done.*");
		if(!($stderr[0] and $stdout[0] and $log[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' ncRNA.$sample.$chk.sh.e* >> $prev/summary.error.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "$prev/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$srcfolder/$chk/$chk.fastq") and (!(-s "$srcfolder/$chk/$chk.R1.fastq") or !(-s "$srcfolder/$chk/$chk.R2.fastq"))){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					#system("rm -f ncRNA.$chk.*");
					#system("rm -f $prev/$chk.done.log");
					system("echo 'Resubmitting the job: ncRNA.$chk.sh' >> job.monitor.txt");
					#system("qsub 00.script/04.ncRNA.script/ncRNA.$chk.sh");
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

## start running the script
opendir(SRC,"$srcfolder") || die"Cannot open file:$!";
my @sams = sort(grep(/\w+/, readdir(SRC)));
closedir SRC;

system("mv ncRNA.$sample.* $prev/");
system("rm -rf $dir");
system("mkdir -p $dir");

## get seq quality scheme
open(SPL, $special);
my %quality = ();
foreach my $line (<SPL>){
	chomp $line;
	my @lines = split(/\s+/, $line);
	if(not exists $quality{$lines[1]}){
		$quality{$lines[1]} = lc($lines[3]);
	}
}
close SPL;

foreach my $sam(@sams)
{
	my $shell = "$dir/trimmomatic.$sample.$sam.sh";
	open(SHL,">$shell") or die "Cannot write $shell: $!";
	
	my $software = 0;
	if($platform eq "sapelo"){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N trimmomatic.$sample.$sam.sh\n";
		print SHL "#PBS -l nodes=1:ppn=1:AMD\n";
		print SHL "#PBS -l walltime=4:00:00\n";
		print SHL "#PBS -l mem=4gb\n\n";
		print SHL "module load java/jdk1.8.0_20\n\n";
		$software = "/usr/local/apps/trimmomatic/0.33/trimmomatic-0.33.jar";
		print SHL "cd \$PBS_O_WORKDIR\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
		$software = "/usr/local/trimmomatic/0.32/trimmomatic-0.32.jar";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	
	my $q = "-phred64";
	if(exists $quality{$sam}){
		if($quality{$sam} eq 'sanger' or $quality{$sam} eq 'illumina1.8'){
			$q = "-phred33";
		}
	}else{
		die "Error: quality scheme has not been identified for $sam\n";
	}

	print SHL "mkdir -p $tgtfolder/$sam\n";
	
	if($mode eq "single-end"){
		print SHL "time java -jar $software SE -threads 1 $q $srcfolder/$sam/$sam.fastq $tgtfolder/$sam/$sam.fastq ILLUMINACLIP:TruSeq.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20\n";
        print SHL "echo \$((`wc -l < $tgtfolder/$sam/$sam.fastq` / 4 )) > $tgtfolder/$sam/$sam.count.txt\n";
	}elsif($mode eq "paired-end"){
		print SHL "time java -jar $software PE -threads 1 $q $srcfolder/$sam/$sam.R1.fastq $srcfolder/$sam/$sam.R2.fastq $tgtfolder/$sam/$sam.R1.fastq $tgtfolder/$sam/$sam.unpaired.R1.fastq $tgtfolder/$sam/$sam.R2.fastq $tgtfolder/$sam/$sam.unpaired.R2.fastq ILLUMINACLIP:TruSeq.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20\n";
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

system("echo 'Finished 05.trimmomatic.pl!' >> job.monitor.txt");
