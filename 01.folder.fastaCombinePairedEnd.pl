#!/usr/bin/perl -w
# run the script: time perl 00.script/01.folder.fastaCombinePairedEnd.pl  01.data/01.Fastq Sapelo

use strict;
system("echo 'Running 01.folder.fastaCombinePairedEnd.pl ....' >> job.monitor.txt");

## read in parameters required by the script
my $srcfolder = shift @ARGV;
my $separater = shift @ARGV || " ";
my $platform = lc(shift @ARGV);
my $sleeptime = shift @ARGV;
my $thread = 1;

my @temp = split(/\//, $srcfolder);
my $sample = pop @temp;
my $prev = "00.script/05.trimmomatic.script/$sample";
my $dir = "00.script/01.checkpairend.script/$sample";
=cut
## check if previous step has succesfully finished
opendir(CHK, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @chks = sort(grep(/^[A-Z0-9]\w+/, readdir(CHK)));
while(1){
        my $count = 0;
        my @temp = @chks;
        my $i = 0;
        while(my $chk = shift @temp){
                my @stderr = glob("trimmomatic.$sample.$chk.sh.o*");
                my @stdout = glob("trimmomatic.$sample.$chk.sh.e*");
                my @log = glob("$prev/$chk.done.*");
                if(!($stderr[0] and $stdout[0] and $log[0])){
                        last; # job hasn't finished, no log file
                }else{
                        system("grep -E 'ERROR|Error|error' trimmomatic.$sample.$chk.sh.e* >> 00.script/05.trimmomatic.script/$sample/summary.error.log\n");
                        $count ++; # job has finished, add one count
                }
        }
        @temp = @chks;
        if($count == scalar @chks){ # all jobs have finished
                if(!-s "00.script/05.trimmomatic.script/$sample/summary.error.log"){ # check if all jobs are successful
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
        }
        sleep $sleeptime;
}
close CHK;
=cut
## start running the script
system("mv trimmomatic.$sample.* $prev/");
system("rm -rf $dir");
system("mkdir -p $dir");

opendir(SRC, $srcfolder) or die "ERROR: Cannot open $srcfolder: $!";
my @subs = sort(grep(/^[A-Z0-9]\w+/, readdir(SRC)));
closedir SRC;
			
foreach my $sub (@subs){
	my $shell = "$dir/fastaCombinePairedEnd.$sample.$sub.sh";
	open(SHL, ">$shell") or die "ERROR: Cannot write $shell: $!";
	
	if($platform eq "sapelo"){
	    print SHL "#PBS -S /bin/bash\n";
	    print SHL "#PBS -q batch\n";
	    print SHL "#PBS -N fastaCombinePairedEnd.$sample.$sub.sh\n";
	    print SHL "#PBS -l nodes=1:ppn=$thread:AMD\n";
	    print SHL "#PBS -l walltime=6:00:00\n";
	    print SHL "#PBS -l mem=20gb\n";
		print SHL "\n";
		print SHL "cd \$PBS_O_WORKDIR\n";
	}elsif($platform eq "zcluster"){
		print SHL "#!/bin/bash\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

	if($platform eq "sapelo"){
		print SHL "module load anaconda/2.2.0\n";
	}elsif($platform eq "zcluster"){
		print SHL "export PYTHONPATH=/usr/local/python/2.7.8/lib/python2.7:/usr/local/python/2.7.8/lib/python2.7/site-packages\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}
    print SHL "time python2.7 00.script/01.fastaCombinePairedEnd.py $srcfolder/$sub/$sub.R1.fastq $srcfolder/$sub/$sub.R2.fastq $separater\n";
    #print SHL "mv $srcfolder/$sub/$sub.R1.fastq_pairs_R1.fastq $srcfolder/$sub/$sub.R1.fastq\n";
    #print SHL "mv $srcfolder/$sub/$sub.R2.fastq_pairs_R2.fastq $srcfolder/$sub/$sub.R2.fastq\n";
    #print SHL "touch $dir/$sub.done.log\n";
	
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
close SRC;

system("echo 'Finished 01.folder.fastaCombinePairedEnd.pl!' >> job.monitor.txt");

