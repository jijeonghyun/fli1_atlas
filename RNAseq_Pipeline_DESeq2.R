#Bulk RNAseq Essentials
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("tidyverse")
BiocManager::install("pheatmap")
BiocManager::install("ggrepel")
BiocManager::install("edgeR")

library("DESeq2")
library("tidyverse")
library("pheatmap")
library("ggrepel")
library("edgeR")

#PLEASE READ BEFORE YOU BEGIN#
#Go to Edit -> Folding -> Collapse All. Go through this code bit by bit.
#Any time you see "Control" and "Treated", you must REPLACE these with the actual sample conditions you have.

#Making a DESeq object----
#read in the raw counts file.
data <- read.csv("raw_counts.csv", row.names = "ensembl_id")
#sort the column names by sample. combine_STAR_output.py will also do this for you.
data <- data[,sort(colnames(data))]

# identify biological replicates. Input your conditions and number of replicates. The words you put in here must match your sample IDs.
condition <- c(rep("Mef2c_KO", 3), rep("NTC", 3))


# assign replicates to each sample name to construct colData
my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(data)
#Check to see if your conditions and sample IDs match up.
my_colData

#Make your DESeq object
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = my_colData,
                              design = ~condition)
dds <- DESeq(dds)

#Normalize and annotate the counts with gene names, then export as a csv file.
normalized_counts <- counts(dds, normalized = T)
annotation <- read.csv("GRCm39.p13_annotation.csv", header = T, stringsAsFactors = F)
normalized_counts <- rownames_to_column(as.data.frame(normalized_counts), var = "ensembl_id")
annotated_data <- right_join(annotation, normalized_counts, by = c("Gene.stable.ID" = "ensembl_id"))
write.csv(annotated_data, file = "gene_annotated_normalized_counts.csv")

#Sample variability analysis----
#Perform a Variance Stabilizing Transformation to analyze sample-to-sample variability
vsd <- vst(dds, blind = TRUE)

#Euclidian distance plot
plotDists = function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste( vsd.obj$condition )
  colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows = sampleDists,
                     clustering_distance_cols = sampleDists,
                     col = colors)
}
plotDists(vsd)

#Variable gene heatmap (NOTE: this is not your DEG heatmap.)
variable_gene_heatmap <- function (vsd.obj, num_genes = 500, annotation, title = "") {
  brewer_palette <- "RdBu"
  # Ramp the color in order to get the scale.
  ramp <- colorRampPalette( RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]
  # get the stabilized counts from the vsd object
  stabilized_counts <- assay(vsd.obj)
  # calculate the variances by row(gene) to find out which genes are the most variable across the samples.
  row_variances <- rowVars(stabilized_counts)
  # get the top most variable genes
  top_variable_genes <- stabilized_counts[order(row_variances, decreasing=T)[1:num_genes],]
  # subtract out the means from each row, leaving the variances for each gene
  top_variable_genes <- top_variable_genes - rowMeans(top_variable_genes, na.rm=T)
  # replace the ensembl ids with the gene names
  gene_names <- annotation$Gene.name[match(rownames(top_variable_genes), annotation$Gene.stable.ID)]
  rownames(top_variable_genes) <- gene_names
  # reconstruct colData without sizeFactors for heatmap labeling
  coldata <- as.data.frame(vsd.obj@colData)
  coldata$sizeFactor <- NULL
  # draw heatmap using pheatmap
  pheatmap::pheatmap(top_variable_genes, color = mr, annotation_col = coldata, fontsize_col = 8, fontsize_row = 250/num_genes, border_color = NA, main = title)
}
variable_gene_heatmap(vsd, num_genes = 50, annotation = annotation)

#Make a PCA plot
plot_PCA = function (vsd.obj) {
  pcaData <- plotPCA(vsd.obj,  intgroup = c("condition"), returnData = T)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    labs(x = paste0("PC1: ",percentVar[1],"% variance"),
         y = paste0("PC2: ",percentVar[2],"% variance"),
         title = "PCA Plot colored by condition") +
    ggrepel::geom_text_repel(aes(label = name), color = "black")
}
plot_PCA(vsd)

#Make a DESeq object comparing two specific conditions (optional)----
#This is only useful if you have >2 sample conditions. If you just have two conditions in your analysis, you can skip this part.
generate_DESeq_object <- function (my_data, groups) {
  data_subset1 <- my_data[,grep(str_c("^", groups[1]), colnames(my_data))]
  data_subset2 <- my_data[,grep(str_c("^", groups[2]), colnames(my_data))]
  my_countData <- cbind(data_subset1, data_subset2)
  condition <- c(rep(groups[1],ncol(data_subset1)), rep(groups[2],ncol(data_subset2)))
  my_colData <- as.data.frame(condition)
  rownames(my_colData) <- colnames(my_countData)
  print(my_colData)
  dds <- DESeqDataSetFromMatrix(countData = my_countData,
                                colData = my_colData,
                                design = ~ condition)
  dds <- DESeq(dds, quiet = T)
  return(dds)
}
ddsObject <- generate_DESeq_object(data, c("Control", "Treated"))
results(ddsObject, contrast = c("condition", "Control", "Treated"))


#Generate DEG results files----
generate_DE_results <- function (dds, comparisons, pvaluecutoff = 0.05, log2cutoff = 0.5, cpmcutoff = 2) {
  # generate average counts per million metric from raw count data 
  raw_counts <- counts(dds, normalized = F)
  cpms <- enframe(rowMeans(edgeR::cpm(raw_counts)))
  colnames(cpms) <- c("ensembl_id", "avg_cpm")
  
  # extract DESeq results between the comparisons indicated
  res <- results(dds, contrast = c("condition", comparisons[1], comparisons[2]))[,-c(3,4)]
  
  # annotate the data with gene name and average counts per million value
  res <- as_tibble(res, rownames = "ensembl_id")
  # read in the annotation and append it to the data
  my_annotation <- read.csv("GRCh38.p13_annotation.csv", header = T, stringsAsFactors = F)
  res <- left_join(res, my_annotation, by = c("ensembl_id" = "Gene.stable.ID"))
  # append the average cpm value to the results data
  res <- left_join(res, cpms, by = c("ensembl_id" = "ensembl_id"))
  
  # combine normalized counts with entire DE list
  normalized_counts <- round(counts(dds, normalized = TRUE),3)
  pattern <- str_c(comparisons[1], "|", comparisons[2])
  combined_data <- as_tibble(cbind(res, normalized_counts))
  
  #Not sure what this does but it was in the original code.
  #combined_data <- as_tibble(cbind(res, normalized_counts[,grep(pattern, colnames(normalized_counts))] ))
  
  #Just get all protein-coding genes
  combined_protein_data <- combined_data[which(combined_data$Gene.type == "protein_coding"),]
  write.csv (combined_protein_data, file = paste0(comparisons[1], "_vs_", comparisons[2], "_all_protein_genes.csv"), row.names =F)
  
  # make ordered rank file for GSEA, selecting only protein coding genes
  res_prot <- res[which(res$Gene.type == "protein_coding"),]
  res_prot_ranked <- res_prot[order(res_prot$log2FoldChange, decreasing = T),c("Gene.name", "log2FoldChange")]
  res_prot_ranked <- na.omit(res_prot_ranked)
  res_prot_ranked$Gene.name <- str_to_upper(res_prot_ranked$Gene.name)
  
  # generate sorted lists with the indicated cutoff values
  res <- res[order(res$log2FoldChange, decreasing=TRUE ),]
  de_genes_pvalue <- res[which(res$pvalue < pvaluecutoff),]
  de_genes_log2f <- res[which(abs(res$log2FoldChange) > log2cutoff & res$padj < pvaluecutoff),]
  de_genes_cpm <- res[which(res$avg_cpm > cpmcutoff & res$padj < pvaluecutoff),]
  
  # write output to files
  write.csv (de_genes_pvalue, file = paste0(comparisons[1], "_vs_", comparisons[2], "_pvalue_cutoff.csv"), row.names =F)
  write.csv (de_genes_log2f, file = paste0(comparisons[1], "_vs_", comparisons[2], "_log2f_cutoff.csv"), row.names =F)
  write.csv (de_genes_cpm, file = paste0(comparisons[1], "_vs_", comparisons[2], "_cpm_cutoff.csv"), row.names =F)
  write.csv (combined_data, file = paste0(comparisons[1], "_vs_", comparisons[2], "_allgenes.csv"), row.names =F)
  write.table (res_prot_ranked, file = paste0(comparisons[1], "_vs_", comparisons[2], "_rank.rnk"), sep = "\t", row.names = F, quote = F)
  
  writeLines( paste0("For the comparison: ", comparisons[1], "_vs_", comparisons[2], ", out of ", nrow(combined_data), " genes, there were: \n", 
                     nrow(de_genes_pvalue), " genes below pvalue ", pvaluecutoff, "\n",
                     nrow(de_genes_log2f), " genes below pvalue ", pvaluecutoff, " and above a log2FoldChange of ", log2cutoff, "\n",
                     nrow(de_genes_cpm), " genes below pvalue ", pvaluecutoff, " and above an avg cpm of ", cpmcutoff, "\n",
                     "Gene lists ordered by log2fchange with the cutoffs above have been generated.") )
  gene_count <- tibble (cutoff_parameter = c("pvalue", "log2fc", "avg_cpm" ), 
                        cutoff_value = c(pvaluecutoff, log2cutoff, cpmcutoff), 
                        signif_genes = c(nrow(de_genes_pvalue), nrow(de_genes_log2f), nrow(de_genes_cpm)))
  invisible(gene_count)
}

#You can edit these cutoffs however you want.
#If you made a DESeq object above, replace "dds" with "ddsObject".
DE_output <- generate_DE_results (dds, c("Mef2c_KO", "NTC"), padjcutoff = 0.05, log2cutoff = 0.5, cpmcutoff = 2)
#The function above will make and export a few more csv files. These will contain all the results (log2foldchange, padj, average cpm, normalized counts, etc.) for each gene, as well as some csv files where the cutoffs were applied.


#Generate a DEG Heatmap----
#Go to your working directory and input the name of your unique resultant filename. It should end with "allgenes.csv".
res <- read.csv ("MEF2C_KO_vs_NTC_allgenes.csv", header = T)

DE_gene_heatmap <- function(res, pvalue_cutoff = 0.05, ngenes = 40) {
  # generate the color palette
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]
  # obtain the significant genes and order by log2FoldChange
  significant_genes <- res %>% filter(pvalue < pvalue_cutoff) %>% arrange (desc(log2FoldChange)) %>% head (ngenes)
  heatmap_values <- as.matrix(significant_genes[,-c(1:7)])
  rownames(heatmap_values) <- significant_genes$Gene.name
  # plot the heatmap using pheatmap
  pheatmap::pheatmap(heatmap_values, color = mr, scale = "row", fontsize_col = 10, fontsize_row = 200/ngenes, fontsize = 5, border_color = NA)
}

#Again, you can adjust your padj cutoff and ngenes.
DE_gene_heatmap(res, 0.005, 100)

#Function to plot single DEGs----
#If you're interested at looking at a specific gene, you can use this function to plot it.
plot_counts <- function (dds, gene, normalization = "DESeq2"){
  # read in the annotation file
  annotation <- read.csv("GRCm39.p13_annotation.csv", header = T, stringsAsFactors = F)
  # obtain normalized data
  if (normalization == "cpm") {
    normalized_data <- cpm(counts(dds, normalized = F)) # normalize the raw data by counts per million
  } else if (normalization == "DESeq2")
    normalized_data <- counts(dds, normalized = T) # use DESeq2 normalized counts
  # get sample groups from colData
  condition <- dds@colData$condition
  # get the gene name from the ensembl id
  if (is.numeric(gene)) { # check if an index is supplied or if ensembl_id is supplied
    if (gene%%1==0 )
      ensembl_id <- rownames(normalized_data)[gene]
    else
      stop("Invalid index supplied.")
  } else if (gene %in% annotation$Gene.name){ # check if a gene name is supplied
    ensembl_id <- annotation$Gene.stable.ID[which(annotation$Gene.name == gene)]
  } else if (gene %in% annotation$Gene.stable.ID){
    ensembl_id <- gene
  } else {
    stop("Gene not found. Check spelling.")
  }
  expression <- normalized_data[ensembl_id,]
  gene_name <- annotation$Gene.name[which(annotation$Gene.stable.ID == ensembl_id)]
  
  #Summarize counts as a tibble for statistics and graphing
  gene_tib <- tibble(condition = condition, expression = expression)
  
  # Plot the DEG, show points and standard error
  ggplot(gene_tib, aes(x = condition, y = expression))+
    geom_jitter(aes(color = condition), size = 3, width = 0.05)+
    stat_summary(fun.data = mean_se, geom = "errorbar", width=0.2)+
    labs (title = paste0("Expression of ", gene_name, " - ", ensembl_id), x = "group", y = paste0("Normalized expression (", normalization , ")"))+
    theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11))
}
#Add any gene name or Ensembl ID and this will make a plot. Error bar is mean +/- standard error.
plot_counts(dds, "Egr1")







#Generate a volcano plot----
res <- read.csv ("Control_vs_Treated_allgenes.csv", header = T)
#This function isn't quite working yet, but you can try it.
plot_volcano <- function (res, padj_cutoff, nlabel = 10, label.by = "padj"){
  # assign significance to results based on padj
  res <- mutate(res, significance=ifelse(res$padj<padj_cutoff, paste0("padj < ", padj_cutoff), paste0("padj > ", padj_cutoff)))
  res = res[!is.na(res$significance),]
  significant_genes <- res %>% filter(significance == paste0("padj < ", padj_cutoff))
  
  # get labels for the highest or lowest genes according to either padj or log2FoldChange
  if (label.by == "padj") {
    top_genes <- significant_genes %>% arrange(padj) %>% head(nlabel)
    bottom_genes <- significant_genes %>% filter (log2FoldChange < 0) %>% arrange(padj) %>% head (nlabel)
  } else if (label.by == "log2FoldChange") {
    top_genes <- head(arrange(significant_genes, desc(log2FoldChange)),nlabel)
    bottom_genes <- head(arrange(significant_genes, log2FoldChange),nlabel)
  } else
    stop ("Invalid label.by argument. Choose either padj or log2FoldChange.")
  
  ggplot(res, aes(log2FoldChange, -log(padj))) +
    geom_point(aes(col=significance)) + 
    scale_color_manual(values=c("red", "black")) + 
    ggrepel::geom_text_repel(data=top_genes, aes(label=head(Gene.name,nlabel)), size = 3)+
    ggrepel::geom_text_repel(data=bottom_genes, aes(label=head(Gene.name,nlabel)), color = "#619CFF", size = 3)+
    labs ( x = "Log2FoldChange", y = "-(Log normalized p-value)")+
    geom_vline(xintercept = 0, linetype = "dotted")+
    theme_minimal()
}
#Again, you can adjust the padjcutoff or the number of labeled genes.
plot_volcano(res, 0.05, 10)
