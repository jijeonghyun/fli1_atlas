# Adapted from wherrylab ATAC pipeline: https://github.com/wherrylab/jogiles_ATAC

# PART ONE #

gunzip *.gz

### Miniconda
# Turn on atacseq miniconda environment
source ~/miniconda/bin/activate
conda activate cutntag

### Load needed modules
module load bowtie2
module load bedtools
module load samtools
module load R

cores=8
memory="8G"
ref="/u/home/j/johnnyji/ref_genome/mm10"
spikeInRef="/u/home/j/johnnyji/ref_genome/ecoli"
chromSize="/u/home/j/johnnyji/mm10.chrom.sizes"
seacr="/u/home/j/johnnyji/SEACR/SEACR_1.3.sh"

### FastQC
# Can check files while alignment is going

#####################################################################

# Run FastQC on all files to get QC stats on raw files
# Change -t # for number of cores
#mkdir fastqc
#find -name '*.fastq' | xargs fastqc -t 8 -o fastqc/

#Trimming off Nextera adapters
echo 'Trimming files...'
for i in *_R1.fastq;
do
 sampleID="${i%_*.*}"
 echo "${sampleID}"
 java -jar /u/home/j/johnnyji/trimmomatic/trimmomatic-0.39.jar PE -threads 8 ${sampleID}_R1.fastq ${sampleID}_R2.fastq ${sampleID}_R1.trimmed.fastq ${sampleID}_R1.untrimmed.fastq ${sampleID}_R2.trimmed.fastq ${sampleID}_R2.untrimmed.fastq ILLUMINACLIP:/u/home/j/johnnyji/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:40 LEADING:3 TRAILING:3
done
echo $?
echo 'Trimming complete!'

mkdir untrimmed
mv *untrimmed.fastq untrimmed

mkdir raw_fastqs
mv *_R1.fastq raw_fastqs
mv *_R2.fastq raw_fastqs


### Alignment
# Mapping cutntag paired end files using bowtie2
# Move bowtie index files into working directory
# Correct the number of cores "p"

echo 'Aligning reads with bowtie2...'
for i in *_R1.trimmed.fastq;
do
 sampleID="${i%_*.*.*}"
 echo "${sampleID}"
 bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 ${sampleID}_R1.trimmed.fastq -2 ${sampleID}_R2.trimmed.fastq -S ${sampleID}.sam &> ${sampleID}.out 
done
echo $?
echo 'Reads aligned with bowtie2!'

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

######################################################################

# PART TWO #

######################################################################

### Make bam files
# Change '@' for number of cores
# Change 'm' for amount of memory PER CORE

echo 'Converting sam to mapped bam...'
for i in `ls *.sam`;
do
	sampleID="${i%.*}"
 	echo "${sampleID}"
	samtools view -bS -F 0x04 -@ ${cores} ${sampleID}.sam >${sampleID}_mapped.bam
done
echo $?
echo 'Converting sam to mapped bam complete!'

echo 'Sorting mapped bam...'
for i in `ls *_mapped.bam`;
do
	sampleID="${i%_*.*}"
 	echo "${sampleID}"
	samtools sort -@ ${cores} -m ${memory} -o ${sampleID}_sorted.bam ${sampleID}_mapped.bam
done
echo $?
echo 'Sorting bam files complete!'

# Make flagstat file to check conversion to bam
for i in `ls *_sorted.bam`;
do
	echo $i >> combFlagstat_bam.txt
	samtools flagstat $i >> combFlagstat_bam.txt
done
echo $?

# Make index for sorted bam files
for i in `ls *_sorted.bam`;
do
 samtools index $i
done
echo $?
echo 'Indices for sorted bam files complete!'

# Make sorted bam files into bigwig file for visualization
echo 'Generating raw bigwig files...'
for i in `ls *_sorted.bam`;
do
	sampleID="${i%_*.*}"
 	echo "${sampleID}"
	bamCoverage --binSize 10 --extendReads -p max -b ${sampleID}_sorted.bam -o ${sampleID}_raw.bw
done
echo $?
echo 'Raw bigwig files complete!'

# STOPPING POINT -- Check bigwig files for visualization

######################################################################

# PART THREE #

######################################################################

### File conversion
# bam -> bed -> fragments bed -> bedgraph

echo 'Converting bam to bed...'
# Bam to bedpe
for i in *_mapped.bam;
do
	sampleID="${i%_*.*}"
	echo "${sampleID}"
	bedtools bamtobed -i ${sampleID}_mapped.bam -bedpe >${sampleID}_bedpe.bed
done
echo $?
echo 'Bam to bed complete!'

# Keep read pairs on the same chromosome and frag length less than 1000bp.
echo 'Cleaning up bedfiles...'
for i in *_bedpe.bed;
do
	sampleID="${i%_*.*}"
	echo "${sampleID}"
	awk '$1==$4 && $6-$2 < 1000 {print $0}' ${sampleID}_bedpe.bed >${sampleID}_clean.bed
done

# Extract fragment related columns
for i in *_clean.bed;
do
	sampleID="${i%_*.*}"
	echo "${sampleID}"
	cut -f 1,2,3,6 ${sampleID}_clean.bed | sort -k1,1 -k2,2n -k3,3n  >${sampleID}_fragments.bed
done

# Remove mitochondrial, random, and non-specific regions
for i in *_fragments.bed;
do
	sampleID="${i%_*.*}"
	echo "${sampleID}"
	sed -i '/chrM/d;/random/d;/chrUn/d' ${sampleID}_fragments.bed 
done
echo $?
echo 'Bedfiles cleaned up!'

#Make bedgraph file for peakcalling
echo 'Converting cleaned up bedfiles to bedgraph...'
for i in *_fragments.bed;
do
 sampleID="${i%_*.*}"
 echo "${sampleID}"
 bedtools genomecov -bg -i ${sampleID}_fragments.bed -g $chromSize > ${sampleID}_fragments.bedgraph
done
echo $?
echo 'Conversion to bedgraph complete!'

######################################################################

### Make folders and sort files

mkdir output

mkdir trimmed_fastqs
mkdir sams
mkdir sortBams
mkdir mapBams
mkdir beds
mkdir bedgraphs
mkdir bws


mv *trimmed.fastq trimmed_fastqs
mv *sam sams
mv *out sams
mv *sorted.bam *sorted.bam.bai sortBams
mv *mapped.bam *mapped.bam.bai mapBams
mv *.bw bws
mv *.bed beds
mv *.bedgraph bedgraphs

mv raw_fastqs output
mv untrimmed output
mv trimmed_fastqs output
mv sams output
mv sortBams output
mv mapBams output
mv beds output
mv bedgraphs output
mv bws output



