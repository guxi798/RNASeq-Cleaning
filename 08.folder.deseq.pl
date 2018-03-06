#!/usr/bin/perl -w
#running command: time perl 00.script/08.folder.deseq.pl 01.data/07.cleanRNA/SRA174408 01.data/11.htseq/SRA174408 01.data/00.PriorData/Athaliana_167_TAIR10.gene.gff3 01.data/12.deseq/SRA174408 Sapelo 10

use strict;
system("echo 'Running 08.folder.deseq.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $reffolder = shift @ARGV;
my $srcfolder = shift @ARGV;
my $genome = shift @ARGV;
my $tissuefile = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my $thread = 1;

my $prevstep = "htseq";
my $currstep = "deseq";
my @temp = split(/\//, $srcfolder);
my $sample = pop @temp;
my $prev = "00.script/07.$prevstep.script/$sample";
my $dir = "00.script/08.$currstep.script/$sample";

## check if previous step has succesfully finished
opendir(CHK, $reffolder) or die "ERROR: Cannot open $reffolder: $!";
my @chks = sort(grep(/^[A-Z]\w+/, readdir(CHK)));
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("$prevstep.$sample.$chk.sh.o*");
		my @stdout = glob("$prevstep.$sample.$chk.sh.e*");
		my @log = glob("$prev/$chk.done.*");
		if(!($stderr[0] and $stdout[0] and $log[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' $prevstep.$sample.$chk.sh.e* >> $prev/summary.error.log\n");
			$count ++; # job has finished, add one count
		}
	}
	@temp = @chks;
	if($count == scalar @chks){ # all jobs have finished
		if(!-s "$prev/summary.error.log"){ # check if all jobs are successful
			system("echo 'There is no error for all jobs' >> job.monitor.txt");
			while(my $chk = shift @temp){
				if(!(-s "$srcfolder/$chk/$chk.counts.txt") and (!(-s "$srcfolder/$chk/$chk.R1.counts.txt") or !(-s "$srcfolder/$chk/$chk.R2.counts.txt"))){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					#system("rm -f $prevstep.$chk.*");
					#system("rm -f $prev/$chk.done.log");
					system("echo 'Resubmitting the job: $prevstep.$chk.sh' >> job.monitor.txt");
					#system("qsub 00.script/07.$prevstep.script/$prevstep.$chk.sh");
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
close CHK;

## start running the script
system("mv $prevstep.$sample.* $prev/");
system("rm -rf $dir");
system("mkdir -p $dir");

	my $shell = "$dir/$currstep.$sample.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	if($platform eq "sapelo"){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N $currstep.$sample.sh\n";
		print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
		print SHL "#PBS -l walltime=48:00:00\n";
		print SHL "#PBS -l mem=40gb\n\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
		print SHL "module load R/3.2.3\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	
	if($platform eq "sapelo"){
		print SHL "mkdir -p $tgtfolder\n";
		print SHL "time /usr/local/R/3.2.3/bin/Rscript 00.script/08.deseq2.R $genome $tissuefile $srcfolder $tgtfolder\n";
		
		print SHL "touch $dir/done.log\n";
		close(SHL);
		system("chmod u+x $shell");
		system("qsub $shell");
	}elsif($platform eq "zcluster"){
		print SHL "mkdir -p $tgtfolder\n";
		print SHL "time /usr/local/R/3.2.3/bin/Rscript 00.script/08.deseq2.R $genome $tissuefile $srcfolder $tgtfolder\n";
		
		print SHL "touch $dir/done.log\n";
		close SHL;
		system("chmod u+x $shell");
		system("qsub -q rcc-30d -l mem_total=20g $shell");
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

system("echo 'Finished 08.folder.deseq.pl!' >> job.monitor.txt");

