#!/usr/bin/perl -w
#running command: time perl 00.script/06.folder.tophat.pl 01.data/01.Fastq/SRA174408 01.data/07.cleanRNA 01.data/00.PriorData/Athaliana_167_TAIR10 01.data/09.tophat single-end Sapelo 10

use strict;
system("echo 'Running 06.folder.tophat.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $reffolder = shift @ARGV;
my $srcfolder = shift @ARGV;
my $genome = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $mode = lc(shift @ARGV);
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my $special = shift @ARGV;
my $thread = 4;

my @temp = split(/\//, $srcfolder);
my $sample = pop @temp;
my $prev = "00.script/01.checkpairend.script/$sample";
my $dir = "00.script/06.tophat.script/$sample";
=cut
## check if previous step has succesfully finished
opendir(CHK, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @chks = sort(grep(/^[A-Z0-9]\w+/, readdir(CHK)));
while(1){
	my $count = 0;
	my @temp = @chks;
	my $i = 0;
	while(my $chk = shift @temp){
		my @stderr = glob("fastaCombinePairedEnd.$sample.$chk.sh.o*");
		my @stdout = glob("fastaCombinePairedEnd.$sample.$chk.sh.e*");
		#my @log = glob("$prev/$chk.done.*");
		#if(!($stderr[0] and $stdout[0] and $log[0])){
        if(!($stderr[0] and $stdout[0])){
			last; # job hasn't finished, no log file
		}else{
			system("grep -E 'ERROR|Error|error' fastaCombinePairedEnd.$sample.$chk.sh.e* >> $prev/summary.error.log\n");
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
					#system("rm -f fastaCombinePairedEnd.$chk.*");
					#system("rm -f $prev/$chk.done.log");
					system("echo 'Resubmitting the job: fastaCombinePairedEnd.$chk.sh' >> job.monitor.txt");
					#system("qsub 00.script/05.trimmomatic.script/trimmomatic.$chk.sh");
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
=cut
## start running the script
opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @sams = sort(grep(/^[A-Z0-9].+/, readdir(SRC)));
closedir SRC;

system("mv fastaCombinePairedEnd.$sample.* $prev/");
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

## construct index for reference
=pod
if($platform eq 'sapelo'){
	system("time /usr/local/apps/bowtie2/2.2.4/bin/bowtie2-build -f -q $genome.fa $genome");
}elsif($platform eq 'zcluster'){
	system("time /usr/local/bowtie2/2.2.3/bin/bowtie2-build -f -q $genome.fa $genome");
}
=cut

foreach my $sam (@sams){
	my $shell = "$dir/tophat.$sample.$sam.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	if($platform eq "sapelo"){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N tophat.$sample.$sam.sh\n";
		print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
		print SHL "#PBS -l walltime=8:00:00\n";
		print SHL "#PBS -l mem=10gb\n\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
    	print SHL "module load tophat/2.0.13\n\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n\n";
		print SHL "export PATH=/usr/local/samtools/0.1.19/:\${PATH}\n";
		print SHL "export LD_LIBRARY_PATH=/usr/local/boost/1.54.0/gcc447/lib:/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:\${LD_LIBRARY_PATH}\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	
	my $t = $thread / 2;
	my $q = '--solexa1.3-quals';
	if(exists $quality{$sam}){
		if($quality{$sam} eq 'sanger' or $quality{$sam} eq 'illumina1.8'){
			$q = "--solexa-quals";
		}
	}else{
		die "Error: quality scheme has not been identified for $sam\n";
	}

	my @R = ();
	if($mode eq "single-end"){
		@R = ("$srcfolder/$sam/$sam.fastq");
	}elsif($mode eq "paired-end"){
		@R = ("$srcfolder/$sam/$sam.R1.fastq", "$srcfolder/$sam/$sam.R2.fastq");
	}else{
		die "Error: please specify mode as either 'single-end' or 'paired-end'";
	}
	
	if($platform eq "sapelo"){
		print SHL "mkdir -p $tgtfolder/$sam\n";
		print SHL "time tophat -i 30 -I 10000 -g 1 -F 0.05 $q --min-segment-intron 30 --max-segment-intron 10000 -G $genome.gene.gff3 --transcriptome-index=transcriptome.index -x 1 -p $thread  -o $tgtfolder/$sam $genome ", join(" ", @R), "\n";
		
		print SHL "touch $dir/$sam.done.log\n";
		close(SHL);
		system("chmod u+x $shell");
		#system("qsub $shell");
	}elsif($platform eq "zcluster"){
		print SHL "mkdir -p $tgtfolder/$sam\n";
		print SHL "time /usr/local/tophat/2.0.13/bin/tophat -i 30 -I 10000 -g 1 -F 0.05 $q --min-segment-intron 30 --max-segment-intron 10000 -G $genome.gene.gff3 --transcriptome-index=transcriptome.index -x 1 -p $thread  -o $tgtfolder/$sam $genome ", join(" ", @R), "\n";
		
		print SHL "touch $dir/$sam.done.log\n";
		close SHL;
		system("chmod u+x $shell");
        system("qsub -q rcc-30d -pe thread $t -l mem_total=10g $shell");
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
}

system("echo 'Finished 06.folder.tophat.pl!' >> job.monitor.txt");

