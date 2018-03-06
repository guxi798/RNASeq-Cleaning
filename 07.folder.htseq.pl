#!/usr/bin/perl -w
#running command: time perl 00.script/07.folder.htseq.pl 01.data/07.cleanRNA/SRA174408 01.data/09.tophat/SRA174408 01.data/00.PriorData/Athaliana_167_TAIR10 01.data/11.htseq/SRA174408 Sapelo 10

use strict;
system("echo 'Running 07.folder.htseq.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $reffolder = shift @ARGV;
my $srcfolder = shift @ARGV;
my $genome = shift @ARGV;
my $tgtfolder = shift @ARGV;
my $mode = shift @ARGV;
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my $thread = 1;

my $prevstep = "tophat";
my $currstep = "htseq";
my @temp = split(/\//, $srcfolder);
my $sample = pop @temp;
$sample = "01.Fastq";
my $prev = "00.script/06.tophat.script/$sample";
my $dir = "00.script/07.$currstep.script/$sample";
=cut
## check if previous step has succesfully finished
opendir(CHK, $reffolder) or die "ERROR: Cannot open $srcfolder: $!";
my @chks = sort(grep(/^[A-Z0-9]\w+/, readdir(CHK)));
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
				if(!(-s "$srcfolder/$chk/accepted_hits.bam")){
					system("echo 'There is no output file for $chk' >> job.monitor.txt");
					#system("rm -f $prevstep.$chk.*");
					#system("rm -f $prev/$chk.done.log");
					system("echo 'Resubmitting the job: $prevstep.$chk.sh' >> job.monitor.txt");
					#system("qsub 00.script/06.$prevstep.script/$prevstep.$chk.sh");
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
opendir(QRY, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @sams = sort(grep(/^\w.+/, readdir(QRY)));
closedir QRY;

system("mv $prevstep.$sample.* $prev/");
system("rm -rf $dir");
system("mkdir -p $dir");

foreach my $sam (@sams){
	my $shell = "$dir/$currstep.$sample.$sam.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	if($platform eq "sapelo"){
		print SHL "#PBS -S /bin/bash\n";
		print SHL "#PBS -q batch\n";
		print SHL "#PBS -N $currstep.$sample.$sam.sh\n";
		print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
		print SHL "#PBS -l walltime=48:00:00\n";
		print SHL "#PBS -l mem=40gb\n\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
    	print SHL "module load samtools/1.2\n";
		print SHL "module load python/2.7.8\n";
		print SHL "module load htseq/0.6.1p1\n\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n\n";
		print SHL "export PYTHONPATH=\${PYTHONPATH}:/usr/local/htseq/0.6.1p1/lib/python\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
	
	if($platform eq "sapelo"){
		print SHL "time samtools sort -n $srcfolder/$sam/accepted_hits.bam $srcfolder/$sam/$sam.sort\n";
        if($mode eq "paired-end"){
            print SHL "samtools view -bf 1 $srcfolder/$sam/$sam.sort.bam > $srcfolder/$sam/htseq_hits.bam\n";
        }else{
            print SHL "mv $srcfolder/$sam/$sam.sort.bam $srcfolder/$sam/htseq_hits.bam\n";
        }
        print SHL "samtools flagstat $srcfolder/$sam/$sam.sort.bam\n";
		print SHL "mkdir -p $tgtfolder/$sam\n";
		print SHL "time htseq-count -f bam -r name -s no -a 10 -t gene -i ID -m union $srcfolder/$sam/htseq_hits.bam $genome.gene.gff3 > $tgtfolder/$sam/$sam.counts.txt\n";
		
		print SHL "touch $dir/$sam.done.log\n";
		close(SHL);
		system("chmod u+x $shell");
		system("qsub $shell");
	}elsif($platform eq "zcluster"){
        print SHL "time /usr/local/samtools/1.2/bin/samtools sort -n $srcfolder/$sam/accepted_hits.bam $srcfolder/$sam/$sam.sort\n";
        if($mode eq "paired-end"){
            print SHL "samtools view -bf 1 $srcfolder/$sam/$sam.sort.bam > $srcfolder/$sam/htseq_hits.bam\n";
        }else{
            print SHL "mv $srcfolder/$sam/$sam.sort.bam $srcfolder/$sam/htseq_hits.bam\n";
        }
        print SHL "samtools flagstat $srcfolder/$sam/$sam.sort.bam\n";
		print SHL "mkdir -p $tgtfolder/$sam\n";
		print SHL "time /usr/local/htseq/0.6.1p1/bin/htseq-count -f bam -r name -s no -a 10 -t gene -i ID -m union $srcfolder/$sam/htseq_hits.bam $genome.gene.gff3 > $tgtfolder/$sam/$sam.counts.txt\n";
		
		print SHL "touch $dir/$sam.done.log\n";
		close SHL;
		system("chmod u+x $shell");
		system("qsub -q rcc-30d -pe thread $thread -l mem_total=10g $shell");
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
}

system("echo 'Finished 07.folder.htseq.pl!' >> job.monitor.txt");

