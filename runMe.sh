#PBS -S /bin/bash
#PBS -q batch
#PBS -N populus_rnaseq
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=12:00:00
#PBS -l mem=4gb
cd $PBS_O_WORKDIR

platform="Sapelo"
mode="paired-end"

:> job.monitor.txt

## decide quality scheme
#time perl 00.script/a2.identify.quality.scheme.pl 01.data/00.PriorData/tissue_record0.txt 01.data/01.Fastq 01.data/00.PriorData/tissue_record.txt 1 10

## process RNAseq data
# 00 select ncRNA for target species, this step only need to be run once
#time perl 00.script/00.select.ncRNA.pl 01.data/00.PriorData 01.data/00.PriorData/total.RNA.fasta Arabidopsis

# 01 cut adapter
time perl 00.script/01.folder.CutAdapter.pl 01.data/01.Fastq 01.data/01.Fastq 01.data/03.CutAdapter $mode $platform 10

# 02 convert fastq to fas
time perl 00.script/02.folder.fas.pl 01.data/01.Fastq 01.data/03.CutAdapter 01.data/04.Fas $mode $platform 10

# 03 blat
time perl 00.script/03.folder.blat.pl 01.data/01.Fastq 01.data/04.Fas 01.data/05.Blat $mode $platform 10

# 04 remove ncRNA
time perl 00.script/04.folder.ncRNA.pl 01.data/01.Fastq 01.data/03.CutAdapter 01.data/05.Blat 01.data/06.ncRNA $mode $platform 10

# 05 remove low-quality reads
time perl 00.script/05.trimmomatic.pl 01.data/01.Fastq 01.data/06.ncRNA 01.data/07.cleanRNA $mode $platform 10 01.data/00.PriorData/tissue_record.txt
time perl 00.script/a3.read.number.summary.pl 01.data/01.Fastq 01.data/03.CutAdapter 01.data/06.ncRNA 01.data/07.cleanRNA $mode $platform 10
time perl 00.script/01.folder.fastaCombinePairedEnd.pl 01.data/07.cleanRNA/$sample " " $platform 10

# 06 align reads to reference
time perl 00.script/06.folder.tophat.pl 01.data/01.Fastq 01.data/07.cleanRNA 01.data/00.PriorData/Pta717s_v1.1 01.data/09.tophat $mode $platform 10 01.data/00.PriorData/tissue_record.txt

# 07 HTseq
time perl 00.script/07.folder.htseq.pl 01.data/01.Fastq 01.data/09.tophat 01.data/00.PriorData/Pta717s_v1.1 01.data/11.htseq $platform 10

# 08 DEseq
time perl 00.script/08.folder.deseq.pl 01.data/01.Fastq 01.data/11.htseq 01.data/00.PriorData/Athaliana_167_TAIR10.gene.gff3 01.data/00.PriorData 01.data/12.deseq $platform 10


# 13 bowtie to map to AS file
#module load bowtie2/2.2.4
#time /usr/local/bowtie2/2.2.3/bin/bowtie2-build -f -q 01.data/00.PriorData/PhQ.transcript.fa 01.data/00.PriorData/PhQ.transcript.fa
#time perl 00.script/13.folder.bowtie.pl 01.data/00.PriorData 01.data/07.cleanRNA 01.data/00.PriorData/PhQ.transcript.fa 01.data/13.bowtie $platform $mode

# retrieve reads mapped to AS junction
#time perl 00.script/a1.folder.extract.read.pl 01.data/13.bowtie 01.data/14.AS.count 01.data/00.PriorData/PhQ.transcript.fa

# 13 abundance estimation using RSEM
#time /usr/local/apps/trinity/r20140717/util/align_and_estimate_abundance.pl --transcripts 01.data/00.PriorData/transcriptome.isoform.fa --est_method RSEM --aln_method bowtie2 --prep_reference
#time perl 00.script/13.estimate.abundance.pl 01.data/07.cleanRNA 13.abundance/whole 01.data/00.PriorData/transcriptome.isoform.fa $mode

#module load R/3.1.2
#Rscript 00.script/a8.summarize.expression.R 13.abundance/whole 13.abundance/whole/gene.fpkm.txt 01.data/00.PriorData/contig2gene.txt 01.data/00.PriorData/tissue_record.txt

