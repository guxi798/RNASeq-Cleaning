#!/usr/bin/perl -w
# run the script: time perl 00.script/01.folder.CutAdapter.pl 01.data/01.Fastq/SRA174408 01.data/01.Fastq/SRA174408 01.data/03.CutAdapter Sapelo 10

use strict;
system("echo 'Running 01.folder.CutAdapter.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $reffolder = shift @ARGV;
my $srcfolder = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $mode = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;

my @temp = split(/\//, $reffolder);
my $sample = pop @temp;
my $prev = "00.script/00.download.script/$sample";
my $dir = "00.script/01.cutadapter.script.count/$sample";

## check if previous step has succesfully finished
=cut
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^\w+/, readdir(CHK)));
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	my @stderr = glob("download.$sample.sh.o*");
	my @stdout = glob("download.$sample.sh.e*");
	my @log = glob("$prev/done.*");
	if(!($stderr[0] and $stdout[0] and $log[0])){
		next; # job hasn't finished, no log file
	}else{
		system("grep -E 'ERROR|Error|error' download.$sample.sh.e* >> $prev/summary.error.log\n");
	}
		if(!-s "$prev/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$srcfolder/$chk/$chk.fastq") and (!(-s "$srcfolder/$chk/$chk.R1.fastq") or !(-s "$srcfolder/$chk/$chk.R2.fastq"))){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
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
	sleep $sleeptime;
}
closedir CHK;
=cut

## start running the script
opendir(SRC,"$srcfolder") || die"Cannot open file:$!";
my @sams = sort(grep(/\w+/, readdir(SRC)));
closedir SRC;

#system("mv download.$sample.* $prev/");
system("rm -rf $dir");
system("mkdir -p $dir");

foreach my $sam(@sams)
{
	my $shell = "$dir/cutadapter.$sample.$sam.sh";
	open(SHL,">$shell") or die "Cannot write $shell: $!";
	
	my $software = 0;
	if($platform eq "sapelo"){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N cutadapter.$sample.$sam.sh\n";
		print SHL "#PBS -l nodes=1:ppn=1:AMD\n";
		print SHL "#PBS -l walltime=12:00:00\n";
		print SHL "#PBS -l mem=4gb\n\n";
		print SHL "module load python/2.7.8\n\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
		$software = "/usr/local/apps/cutadapt/1.9.dev1/bin/cutadapt";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
		$software = "/usr/local/cutadapt/1.9.dev1/bin/cutadapt";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	
	print SHL "mkdir -p $tgtfolder/$sam\n";
	
	if($mode eq "single-end"){
		print SHL "time python2.7 $software -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -O 4 -m 10 -o $tgtfolder/$sam/$sam.fastq $srcfolder/$sam/$sam.fastq  \n";
        print SHL "echo \$((`wc -l < $srcfolder/$sam/$sam.fastq` / 4 )) > $srcfolder/$sam/$sam.count.txt\n";
        print SHL "echo \$((`wc -l < $tgtfolder/$sam/$sam.fastq` / 4 )) > $tgtfolder/$sam/$sam.count.txt\n";
	}elsif($mode eq "paired-end"){
		#print SHL "time python2.7 $software -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -O 4 -m 10  -o $tgtfolder/$sam/$sam.R1.fastq $srcfolder/$sam/$sam.R1.fastq\n";
		#print SHL "time python2.7 $software -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -O 4 -m 10  -o $tgtfolder/$sam/$sam.R2.fastq $srcfolder/$sam/$sam.R2.fastq\n";
        print SHL "echo \$((`wc -l < $srcfolder/$sam/$sam.R1.fastq` / 4 )) > $srcfolder/$sam/$sam.R1.count.txt\n";
        print SHL "echo \$((`wc -l < $srcfolder/$sam/$sam.R2.fastq` / 4 )) > $srcfolder/$sam/$sam.R2.count.txt\n";
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

system("echo 'Finished 01.folder.CutAdapter.pl!' >> job.monitor.txt");
