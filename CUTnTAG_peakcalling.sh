
# Turn on atacseq miniconda environment
source ~/miniconda/bin/activate
conda activate cutntag

### Load needed modules
module load bowtie2
module load bedtools
module load samtools
module load R

cores=8
ref="/u/home/j/johnnyji/ref_genome/mm10"
spikeInRef="/u/home/j/johnnyji/ref_genome/ecoli"
chromSize="/u/home/j/johnnyji/mm10.chrom.sizes"
seacr="/u/home/j/johnnyji/SEACR/SEACR_1.3.sh"


for i in *_fragments.bedgraph;
do
        sampleID="${i%_*.*}"
        echo "${sampleID}"
        bash $seacr ${sampleID}_fragments.bedgraph 0.01 norm stringent ${sampleID}_seacr_norm_top0.01.peaks
done

#I call normalized 0.01% peaks based on AUC.


