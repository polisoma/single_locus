# ES open regions: files are strange and do not look like peaks! --> will not use them
# extend all the peaks to 1000 bp (+/- 500 bp)
# files: 
#ZHBTc4-mESC/ENCFF001YVO	doxycycline hyclate
#ZHBTc4-mESC/ENCFF001YVP	doxycycline hyclate
#ZHBTc4-mESC/ENCFF111WVB	doxycycline hyclate
#ZHBTc4-mESC/ENCFF082FXK	doxycycline hyclate
#ZHBTc4-mESC/ENCFF001YVL	doxycycline hyclate
#ZHBTc4-mESC/ENCFF001YMF	doxycycline hyclate
#ZHBTc4-mESC/ENCFF001YVK	doxycycline hyclate

# ZHBTc4-mESC/ENCFF001YVV.bed.gz
# ZHBTc4-mESC/ENCFF001YMH.bed.gz
# R1/ENCFF267ODA.bed.gz
# ES-CJ7/ENCFF001YKZ.bed.gz  
# ES-CJ7/ENCFF001YOG.bed.gz  
# ES-CJ7/ENCFF001YOH.bed.gz
# ES-E14/ENCFF001YOK.bed.gz  
# ES-E14/ENCFF001YOL.bed.gz  
# ES-E14/ENCFF167GZG.bed.gz
# ES-E14/ENCFF911LAE.bed.gz
##########################################
# files downloaded directly from UCSC browser, HotSpots
# extending all files and merging   
# bedtools: /usr/local/src/BEDTools-Version-2.16.2/bin/mergeBed -h

cd ~/data/published_data/ENCODE/mouse/HotSpots
rm /tmp/all_merged.txt.gz
for file in wgEncodeUwDnaseEscj7S129ME0HotspotsRep1.broadPeak.gz wgEncodeUwDnaseEse14129olaME0HotspotsRep2.broadPeak.gz wgEncodeUwDnaseEscj7S129ME0HotspotsRep2.broadPeak.gz wgEncodeUwDnaseEsww6UknME0HotspotsRep1.broadPeak.gz wgEncodeUwDnaseEse14129olaME0HotspotsRep1.broadPeak.gz wgEncodeUwDnaseEsww6UknME0HotspotsRep2.broadPeak.gz; do
zcat $file | awk '{print $1,$2,$3}' | gzip >> /tmp/all_merged.txt.gz
done

zcat /tmp/all_merged.txt.gz | sort -k1,1 -k2,2n | awk -vOFS="\t" '{print $1,$2,$3}' > /tmp/all_dasha.txt
/usr/local/src/BEDTools-Version-2.16.2/bin/mergeBed -i /tmp/all_dasha.txt > all_merged_final.txt
gzip all_merged_final.txt
rm /tmp/all_dasha.txt /tmp/all_merged.txt.gz

# Making bowtie index for open regions 
# get fasta
bedtools getfasta -fi /home/daria/data/data_genomes/mm9/fa/mm9.fa -s -bed <(zcat all_merged_final.txt.gz) -fo all_merged_final.fa
# making index
bowtie-build -f all_merged_final.fa all_merged_final 

###########################################
cd ~/data/published_data/ENCODE/mouse/HotSpots

# picking sgRNA and mapping them to the "open" genome of mESC
# removing all spaces in a text file
cat DE_sequence.fa | tr -d ' ' | awk '{print ">DE_sequence_Oct4""\n"$0}' > DE_seq.fa

# SSC-based selection of guides (whole regions but )
~/utils/SSC0.1/SSC0.1/bin/Fasta2Spacer -5 20 -3 0 -i sgRNA/DE_seq.fa -o sgRNA/DE_sequence.spacer.seq
~/utils/SSC0.1/SSC0.1/bin/Fasta2Spacer -5 20 -3 0 -i sgRNA/PE_seq.fa -o sgRNA/PE_sequence.spacer.seq

# SSC score
~/utils/SSC0.1/SSC0.1/bin/SSC -l 20 -m ~/utils/SSC0.1/SSC0.1/matrix/human_CRISPRi_20bp.matrix -i sgRNA/DE_sequence.spacer.seq -o sgRNA/DE_sequence.spacer.out
~/utils/SSC0.1/SSC0.1/bin/SSC -l 20 -m ~/utils/SSC0.1/SSC0.1/matrix/human_CRISPRi_20bp.matrix -i sgRNA/PE_sequence.spacer.seq -o sgRNA/PE_sequence.spacer.out

# SSC score in the paper: median ~0.08 => arbitrary cut off
# creat list of sgRNA with all possible PAMs
sort -k6,6gr sgRNA/DE_sequence.spacer.out | awk '($6>=0.05 && ($3<=296 || $2>=452)){print ">seq"$6"\n"$1"GGG""\n"">seq"$6"\n"$1"AGG""\n"">seq"$6"\n"$1"CGG""\n"">seq"$6"\n"$1"TGG"}' > /tmp/DE_sgRNA.fa
sort -k6,6gr sgRNA/PE_sequence.spacer.out | awk '($6>=0.05 && ($3<=446 || $2>=817)){print ">seq"$6"\n"$1"GGG""\n"">seq"$6"\n"$1"AGG""\n"">seq"$6"\n"$1"CGG""\n"">seq"$6"\n"$1"TGG"}' > /tmp/PE_sgRNA.fa

# mapping sgRNA to bowtie open genome index: -a show all alighnments
cat /tmp/DE_sgRNA.fa | bowtie -p 3 -a -v 3 --fullref --quiet ~/data/published_data/ENCODE/mouse/HotSpots/all_merged_final -f -

#### greping on specific stuff
cat /tmp/PE_sgRNA.fa | bowtie -p 3 -a -v 3 --fullref --quiet ~/data/published_data/ENCODE/mouse/HotSpots/all_merged_final -f - | awk '{count[$1]++}END{for(i in count)print i,count[i]}' | sort -k2,2n

#########
 
# copying all files but not directories
cd /Users/tfg900/Dropbox/workBRIC/PROJECTS/DSH-65122057/FASTQ_Generation_2018-02-13_01_16_03Z-79418741

for f in * ; do
echo -en "$f\n"
scp $f/*.fastq.gz daria@tycho.sund.root.ku.dk:~/data/projects/BRIC_data/single_locus/raw/.
done

# mapping with bowtie

# pulling lanes together (bowtie 1 does not take individual lanes)
 
myTemp=~/data/temp_NOBACKUP 
for f in `ls raw/*_R1_001.fastq.gz`; do   
   name=$(echo -en "${f%%.*}" | awk '{split($1,a,"_S"); split(a[1],s,"DSH-"); print s[2]}')
   label=$(basename $name)
   echo -en "$label\n"
   zcat $f >> $myTemp/$label.merged.fastq
done

# number of unmapped reads: 
myTemp=~/data/temp_NOBACKUP
for name in $myTemp/*.merged.fastq; do
	short=$(basename $name .merged.fastq)
	echo -en "$short\n"
	cat $myTemp/$short.merged.fastq | wc -l
done

# devide numbers by 4! 
#DE2-input-rep1
#87798236
#DE2-IP-rep1
#73370156
#DE2-IP-rep2
#101800408
#NC-input-rep1
#74270944
#NC-IP-rep1
#54802048
#NC-IP-rep2
#78215720
#PE3-input-rep1
#73701124
#PE3-IP-rep1
#54044600

screen -R mapping
# do not have gz files - cat on files 
# will map to mm10
mm10=/k/genomes/mm10/index/bowtie_canonical/mm10 # somehow index in my directory was not working
myTemp=~/data/temp_NOBACKUP

for name in $myTemp/*.merged.fastq; do
	short=$(basename $name .merged.fastq)
	echo -en "$short\n"
	cat $myTemp/$short.merged.fastq | bowtie -p 5 -q -m 1 -v 3 --sam --best --strata --quiet $mm10 - > $myTemp/$short.sam
done

# for the future add samtools rmdup -s bam/${file}.bam bam/${file}.nodup.bam on sorted bam --> so all bam files are no dup
myTemp=~/data/temp_NOBACKUP
for name in $myTemp/*.merged.fastq; do
	sample=$(basename $name .merged.fastq)
	echo -en "$sample\n"
	samtools view -Sb $myTemp/${sample}.sam > $myTemp/${sample}_nonSorted.bam
	samtools sort $myTemp/${sample}_nonSorted.bam > $myTemp/${sample}.bam
	samtools index $myTemp/${sample}.bam
done

# after everything is done
myTemp=~/data/temp_NOBACKUP
for name in $myTemp/*.merged.fastq; do
	sample=$(basename $name .merged.fastq)
	rm $myTemp/${sample}.sam $myTemp/${sample}_nonSorted.bam
done

myTemp=~/data/temp_NOBACKUP
mkdir -p ~/data/projects/BRIC_data/single_locus/bam
finalFolder=~/data/projects/BRIC_data/single_locus/bam
for name in $myTemp/*.merged.fastq; do
	sample=$(basename $name .merged.fastq)
	echo -en "$sample\n"
	mv $myTemp/${sample}.bam  $finalFolder/.
	mv $myTemp/${sample}.bam.bai  $finalFolder/.
done 

cd ~/data/projects/BRIC_data/single_locus
mkdir -p bw
myTemp=~/data/temp_NOBACKUP
for name in $myTemp/*.merged.fastq; do
	sample=$(basename $name .merged.fastq)
	echo -en "$sample\n"
  EXTEND=140
  CHROMSIZE=/home/daria/data/data_genomes/mm10/chrom.sizes
  # Number of reads
  librarySize=$(samtools idxstats bam/${sample}.bam | awk '{total+=$3} END{print total}')
  echo $librarySize
  # Create density file: extend reads, calculate read density at each position and normalize the library size to 1 million reads
  bamToBed -i bam/${sample}.bam | \
  awk -vCHROM=$CHROMSIZE -vEXTEND=$EXTEND -vOFS='\t' 'BEGIN{while((getline<CHROM)>0){chromSize[$1]=$2}}{chrom=$1;start=$2;end=$3;strand=$6;if(strand=="+"){end=start+EXTEND;if(end>chromSize[chrom]){end=chromSize[chrom]}};if(strand=="-"){start=end-EXTEND;if(start<1){start=1}};print chrom,start,end}' | \
  sort -k1,1 -k2,2n | genomeCoverageBed -i stdin -g $CHROMSIZE -d | \
  awk -vOFS='\t' -vSIZE=$librarySize '{print $1,$2,$2+1,$3*1000000/SIZE}' | gzip > $myTemp/${sample}.density.gz	
  # Create WIG file
  gunzip -c $myTemp/${sample}.density.gz | \
  awk -vOFS='\t' '($4!=0){if(!chrom[$1]){print "variableStep chrom="$1;chrom[$1]=1};print $2,$4}' | \
  gzip > $myTemp/${sample}.wig.gz
  # Create BigWig file, takes quite some time and memory-consuming
  wigToBigWig $myTemp/${sample}.wig.gz $CHROMSIZE bw/${sample}.bw
  # Remove intermediate file
  rm $myTemp/${sample}.wig.gz $myTemp/${sample}.density.gz
done

########## making bw with macs2
# filtering dup and making bed files

screen -R tracks
cd ~/data/projects/BRIC_data/single_locus
myTemp=~/data/temp_NOBACKUP
mkdir -p macs2 

for name in $myTemp/*.merged.fastq; do
	sample=$(basename $name .merged.fastq)
	echo -en "$sample\n"
   macs2 filterdup -i bam/${sample}.bam --keep-dup=1 -o macs2/${sample}_filterdup.bed
   macs2 pileup -i macs2/${sample}_filterdup.bed -o macs2/${sample}_filterdup.pileup.bdg --extsize 140
done

# bedGraphtoBigWig
# claims that chrM is bigger and when removed - complained: "enlarged" chrM and saved tp the same folder as data
#CHROMSIZE=chrom.sizes # problems with chromoc=somal sizes
#bedGraphToBigWig macs2/DE2-input-rep1_filterdup.pileup.bdg $CHROMSIZE bw/DE2-input-rep1_filterdup.pileup.bw

CHROMSIZE=chrom.sizes 
for name in $myTemp/*.merged.fastq; do
	sample=$(basename $name .merged.fastq)
	echo -en "$sample\n"
	bedGraphToBigWig macs2/${sample}_filterdup.pileup.bdg $CHROMSIZE bw/${sample}_filterdup.pileup.bw
done

# macs2 tracks were much more faster
# making soft links to web folder

cd ~/web/projects/single_locus
folder=/NextGenSeqData/project-data/daria/projects/BRIC_data/single_locus/bw
for file in $folder/*bw; do
	ln -s $file .
done

# making links for all files: 
for file in ~/web/projects/single_locus/*bw; do
	name=$(basename $file)
	echo -en "https://bricweb.sund.ku.dk/bric-data/daria/projects/single_locus/$name\n"
done

# calling peaks with macs2 and then peakzilla
# differential peak calling using NC controls

cd ~/data/projects/BRIC_data/single_locus/
mkdir -p ~/data/projects/BRIC_data/single_locus/macs2/differential

screen -R tracks
# no second replicate for inputs
# by mistake called the final file "input" --> not input but IP!
for rep in rep2; do
	macs2 callpeak -B -t bam/DE2-IP-${rep}.bam -c bam/DE2-input-rep1.bam -n DE2-input-${rep} --outdir macs2/differential --nomodel --extsize 200
	macs2 callpeak -B -t bam/NC-IP-${rep}.bam -c bam/NC-input-rep1.bam -n NC-input-${rep} --outdir macs2/differential --nomodel --extsize 200
done

# differential peak calling
# need to know number of tags after filtering

tags1=11681290 # for DE2
tags2=8730593 # for NC
folder=macs2/differential
# -g and -l parameters left unchanged
macs2 bdgdiff --t1 $folder/DE2-input-rep1_treat_pileup.bdg --c1 $folder/DE2-input-rep1_control_lambda.bdg \
 --t2 $folder/NC-input-rep1_treat_pileup.bdg --c2 $folder/NC-input-rep1_control_lambda.bdg \
 --d1 $tags1 --d2 $tags2 -g 60 -l 120 --o-prefix $folder/diff_DE2-IP-rep1_vs_NC-IP-rep1

# for second replicates

tags1=16839162 # for DE2
tags2=12642723 # for NC
folder=macs2/differential
# -g and -l parameters left unchanged
macs2 bdgdiff --t1 $folder/DE2-input-rep2_treat_pileup.bdg --c1 $folder/DE2-input-rep2_control_lambda.bdg \
 --t2 $folder/NC-input-rep2_treat_pileup.bdg --c2 $folder/NC-input-rep2_control_lambda.bdg \
 --d1 $tags1 --d2 $tags2 -g 60 -l 120 --o-prefix $folder/diff_DE2-IP-rep2_vs_NC-IP-rep2

# counting number of reads for each bam file on a union of peaks from 2 replicates
# take common peaks and unique for D2 each replicate, find summit and take 100 bp only (+/-50 bp)

cat macs2/differential/diff_DE2-IP-rep1_vs_NC-IP-rep1_c3.0_cond1.bed \
macs2/differential/diff_DE2-IP-rep2_vs_NC-IP-rep2_c3.0_cond1.bed | \
awk -vOFS="\t" -vWIN=200 '($1!="track"){summit=$2+int(($3-$2)/2);print $1,summit-WIN,summit+WIN,$4}' | sort -k1,1 -k2,2n > /tmp/allPeaks.bed

# chr17   35504095        35504321

# use bams with filtered duplikates --> bed --> extend the tags --> overlap with a set of peaks --> DESeq
# set of peaks: unique to each of DE2 replictes + common? or a set of random regions?
# in the paper: predicted offtrgets...


####################### NC dcas9 cell line --> non-specific chipable regions with flag

# call peaks on 2 replicates or differential
# overlap or concentrate on common peaks between 2 reps
# how many are there?

# calling peaks with macs2 and then peakzilla
# differential peak calling using NC controls

cd ~/data/projects/BRIC_data/single_locus/
mkdir -p ~/data/projects/BRIC_data/single_locus/macs2/differential_NC_only

screen -R tracks
# no second replicate for inputs
# by mistake called the final file "input" --> not input but IP!
for rep in rep2; do
	macs2 callpeak -B -t bam/DE2-IP-${rep}.bam -c bam/DE2-input-rep1.bam -n DE2-input-${rep} --outdir macs2/differential --nomodel --extsize 200
	macs2 callpeak -B -t bam/NC-IP-${rep}.bam -c bam/NC-input-rep1.bam -n NC-input-${rep} --outdir macs2/differential --nomodel --extsize 200
done

# differential peak calling
# need to know number of tags after filtering

tags1=11681290 # for DE2
tags2=8730593 # for NC
folder=macs2/differential
# -g and -l parameters left unchanged
macs2 bdgdiff --t1 $folder/DE2-input-rep1_treat_pileup.bdg --c1 $folder/DE2-input-rep1_control_lambda.bdg \
 --t2 $folder/NC-input-rep1_treat_pileup.bdg --c2 $folder/NC-input-rep1_control_lambda.bdg \
 --d1 $tags1 --d2 $tags2 -g 60 -l 120 --o-prefix $folder/diff_DE2-IP-rep1_vs_NC-IP-rep1













