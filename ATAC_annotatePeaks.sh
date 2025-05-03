source ~/miniconda/bin/activate
conda activate cutntag

cores=8
ref="/u/home/j/johnnyji/ref_genome/mm10"
chromSize="/u/home/j/johnnyji/mm10.chrom.sizes"

#Merge all peaks together to make consensus peak file
#Depending on what you need, you will want to either --intersect the peaks (takes only what is common across samples)
#Or you will want to --merge the peaks (merges all peaks together)
bedops --intersect *.narrowPeak > masterPeaks.bed

#Convert consensus peak file into homer-readable format
awk '{print $1"\t"$2"\t"$3"\t"NR"\t"$1":"$2"-"$3"\t""+"}' masterPeaks.bed | awk 'BEGIN {print "Chr""\t""Start""\t""End""\t""Peakid""\t""Region""\t""Strand"} {print}' > homerPeaks.bed

#Use homer to annotate the master peak file
#Must change genome if necessary.
annotatePeaks.pl homerPeaks.bed mm10 > annoPeaks.txt

#Convert annotated peak file into saf format
awk 'BEGIN {OFS="\t"} NR>1 {print $1, $2, $3, $4, "."}' annoPeaks.txt > annoPeaks.saf

#Must use -p if you want to specify that it is paired-end.
featureCounts -a annoPeaks.saf -F SAF -o featureCounts.txt *_debl.bam















