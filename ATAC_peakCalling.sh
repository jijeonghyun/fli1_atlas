### Miniconda
# Turn on atacseq miniconda environment
source ~/miniconda/bin/activate
conda activate atacseq

### Load needed modules
module load bowtie2
module load bedtools
module load samtools
module load R

### Alignment
# Mapping atac paired end files to hg38 using bowtie2
# Move bowtie index files into working directory
# Correct the number of cores "p"
cores=8
ref="/u/home/j/johnnyji/ref_genome/mm10"
chromSize="/u/home/j/johnnyji/mm10"

echo 'Calling peaks with MACS2...'
for i in `ls *_debl.bed`;
do
    sampleID="${i%.*}"
    echo "${sampleID}"
    macs2 callpeak -t $i -g mm -n ${sampleID} -q 0.01
done

echo $?
echo 'You have called peaks!'














