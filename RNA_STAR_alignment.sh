# From wherrylab ATAC pipeline: https://github.com/wherrylab/jogiles_ATAC

# PART ONE #

### Miniconda
# Turn on atacseq miniconda environment
source ~/miniconda/bin/activate
conda activate rnaseq


### Alignment
# Mapping rna paired end files to hg38 using STAR
# Move STAR index files into working directory
# Correct the number of cores "p"

mkdir star_output
echo 'Aligning reads with STAR...'
for i in *R1.trimmed.fastq;
do
 sampleID="${i%_*.*.*}"
 echo "${sampleID}"
 STAR --runThreadN 16 --genomeDir $SCRATCH/star_index/mm39 -c --outFilterType BySJout --readFilesIn ${sampleID}_R1.trimmed.fastq ${sampleID}_R2.trimmed.fastq --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix ./star_output/${sampleID}.
done
echo $?
echo 'Alignment complete!'

######################################################################

