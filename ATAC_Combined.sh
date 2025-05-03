# From wherrylab ATAC pipeline: https://github.com/wherrylab/jogiles_ATAC

# PART ONE #

### Miniconda
# Turn on atacseq miniconda environment
source ~/miniconda/bin/activate
conda activate atacseq

### FastQC
# Can check files while alignment is going

######################################################################

# Run FastQC on all files to get QC stats on raw files
# Change -t # for number of cores
#mkdir fastqc
#find -name '*.fastq' | xargs fastqc -t 16 -o fastqc/


cores=8 #This must be just a number
memory="8G" #This must be upper case G
totalram="64g" #this must be lower case g
ref="/u/home/j/johnnyji/ref_genome/mm10"
chromSize="/u/home/j/johnnyji/mm10.chrom.sizes"
blackList="/u/home/j/johnnyji/mm10.blacklist.bed"

### Alignment
# Mapping atac paired end files to mm10 using bowtie2
# Move bowtie index files into working directory
# Correct the number of cores "p"

echo 'Aligning reads using bowtie2...'
for i in *_R1_trimmed.fastq;
do
 sampleID="${i%_*_*.*}"
 echo "${sampleID}"
 bowtie2 -k1 -N1 -p${cores} -x ${ref} -1 ${sampleID}_R1_trimmed.fastq -2 ${sampleID}_R2_trimmed.fastq -S ${sampleID}.sam 2> ${sampleID}.out
done
echo $?
echo 'Read alignment complete!'

# Combine all of the .out files produced from the sam alignment

for i in *.out;
do
	sampleID="${i%.*}"
	echo "${sampleID}" > ${sampleID}_out.txt
	cat $i >> ${sampleID}_out.txt
done

find . -name '*_out.txt' | xargs cat >> combOutSams.txt
rm *_out.txt

# STOPPING POINT -- Check fastqc and alignment reports

echo 'Sorting bam files...'
for i in `ls *.sam`;
do
	sampleID="${i%.*}"
 	echo "${sampleID}"
	samtools view -bS -@ ${cores} $i | samtools sort - -m ${memory} -@ ${cores} -T ${sampleID}_sort -o ${sampleID}_sort.bam
done
echo $?
echo 'Bam files sorted!'


# Make flagstat file to check conversion to bam
for i in `ls *_sort.bam`;
do
	echo $i >> combFlagstat_bam.txt
	samtools flagstat $i >> combFlagstat_bam.txt
done
echo $?


# Make index for sorted bam files
for i in `ls *_sort.bam`;
do
 samtools index $i
done


# Make bam files with only mapped reads
echo 'Extracting mapped-only bam files...'
for i in `ls *_sort.bam`;
do
	sampleID="${i%_*.*}"
 	echo "${sampleID}"
	samtools view -bS -f 2 -@ ${cores} $i > ${sampleID}_map.bam
done
echo $?
echo 'Mapped bam files complete!'

# Make index for mapped-only bam files
for i in `ls *_map.bam`;
do
 samtools index $i
done


# Make mapped bam files without mitochondrial reads
echo 'Removing mitochondrial reads...'
for i in `ls *_map.bam`;
do
    sampleID="${i%*.*}"
    echo "${sampleID}"
    samtools view -b -@ ${cores} $i chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > ${sampleID}_deMito.bam
done
echo $?
echo 'Mitochondrial reads removed!'


# Make index for deMito bams
for i in `ls *_deMito.bam`;
do
 samtools index $i
done

echo 'Sorting bam files again...'
for i in `ls *_deMito.bam`;
do
	sampleID="${i%.*}"
 	echo "${sampleID}"
	samtools sort -m ${memory} -@ ${cores} -T ${sampleID}_sort -o ${sampleID}_sort.bam ${sampleID}.bam
done
echo $?
echo 'Bam files sorted by coordinate!'

######################################################################

### Remove duplicates

# Change 'Xmx_g for amount of available ram

echo 'Removing duplicates...'
for i in `ls *deMito_sort.bam`;
do
	sampleID="${i%_*.*}"
 	echo "${sampleID}"
	picard -Xmx${totalram} MarkDuplicates INPUT=$i OUTPUT=${sampleID}_dedup.bam METRICS_FILE=${sampleID}_dedupMetrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORT_ORDER=coordinate REMOVE_DUPLICATES=true
done
echo $?
echo 'PCR duplicates removed!'


echo 'Making indices...'

# Make index for deduped bam files
for i in `ls *_dedup.bam`;
do
 samtools index $i
done
echo $?
echo 'You have indexes for deduplicated bam files!'

# Combine picard duplicate metrics
for i in `ls *_dedupMetrics.txt`;
do
	sampleID="${i%_*_*_*.*}"
	awk 'FNR==8 {print FILENAME, $0}' ${sampleID}_map_deMito_dedupMetrics.txt >> comb_picardDedupStats1.txt
done

awk '{print $1"\t"$9"\t"$10"\t"$11}' comb_picardDedupStats1.txt > comb_picardDedupStats2.txt


echo $'sampleID\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE' | cat - comb_picardDedupStats2.txt > comb_picardDedupStats.txt
rm comb_picardDedupStats1.txt
rm comb_picardDedupStats2.txt


### Convert bams to beds --> shift reads to reflect center of Tn5 insertion
#(The Tn5 transposon binds as a dimer and inserts two adaptors separated by 9 bp, all reads aligning to the + strand were offset by +4 bp, and all reads aligning to the – strand were offset −5 bp.)

echo 'Converting bam to bed and shifting reads to reflect center of Tn5 insertion...'
for i in `ls *_dedup.bam`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	bedtools bamtobed -i  $i | awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 5, $3 + 5, $4, $5, $6; else print $1, $2 - 4, $3 - 4, $4, $5, $6}' > ${sampleID}.bed
done

echo $?
echo 'Bedfiles complete!'


#####################################################################


### Remove blacklisted regions

# Copy mm10.blacklist.bed into working directory
echo 'Removing blacklist regions from bedfile...'
for i in `ls *.bed`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	bedtools subtract -a $i -b ${blackList} > ${sampleID}_debl.bed
done

echo $?
echo 'Blacklist regions removed!'



#####################################################################


### Make bam files and associated indices from deblacklisted bed files
echo 'Converting bed back to bam...'
for i in `ls *_debl.bed`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	bedtools bedtobam -i $i -g ${chromSize} | samtools sort - -m ${memory} -@ ${cores} -o ${sampleID}.bam
done

echo $?
echo 'Deblackedlisted bam files complete!'

# Make index for deblacklisted bam files
for i in `ls *_debl.bam`;
do
 samtools index $i
done

### Make normalized bigwig files

# Must copy mm10.chrom.size file into working directory
# Turn on wigtobigwig environment in conda

eval "$(conda shell.bash hook)"
conda activate wigtobigwig
echo 'Making bigiwigs...'
for i in `ls *_debl.bam`;
do
 sampleID="${i%.*}"
 lines=$(samtools view -c ${sampleID}.bam);\
 bedtools genomecov -ibam ${sampleID}.bam -bg -scale $(echo "1000000 / ${lines} " | bc -l) -g ${chromSize} | \
 wigToBigWig -clip stdin ${chromSize} ${sampleID}.bw 2> ${sampleID}_log
done

echo $?
echo 'Bigwigs complete!'

# Go back to atacseq environment in conda
conda activate atacseq

#####################################################################

#Make bedgraph file for peakcalling
echo 'Converting cleaned up bedfiles to bedgraph...'
for i in `ls *_debl.bed`;
do
	sampleID="${i%_*.*}"
	echo "${sampleID}"
	bedtools genomecov -bg -i ${sampleID}_debl.bed -g $chromSize > ${sampleID}.bedgraph
done

echo $?
echo 'Conversion to bedgraph complete!'













