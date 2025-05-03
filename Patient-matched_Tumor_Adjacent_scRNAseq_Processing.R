#Load libraries and set the working directory ----------
library(SeuratObject)
library(Seurat)
library(harmony)
library(DESeq2)
library(tidyverse)
library(harmony)
library(ggpubr)
library(scCustomize)

getwd()
setwd("/Users/lab/Documents/areeba/Codes for analysis")
#Function to add new desired columns ----------
addcolumns <- function (dataset, cancertype, celltype, condition, source, PatientID) {
  # Add a cancertype column, then add the cancertype to that column. 
  columns_to_keep <- c("nCount_originalexp", "nFeature_originalexp","PatientID",
                       "cancer_type", "cell_type",
                       "Condition","Source")
  dataset@meta.data$cancer_type <- cancertype
  dataset@meta.data$cell_type <- celltype
  dataset@meta.data$Condition <- condition
  dataset@meta.data$Source <- source
  dataset@meta.data <- dataset@meta.data[, columns_to_keep]
  print(dataset)
  
}
#Cleaning up data ----------
###Running PCA and filtering out T - cells (done on normal and tumor separately)
TEST <- LungN_Clean_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

LungN_Clean_Final <- TEST

DimPlot(LungN_Clean_Final, reduction = "umap", label = TRUE)
VlnPlot(LungN_Clean_Final, features = Tcellgenes)

Idents(Prostate_Low_Grade_Normal_Final) <- 'seurat_clusters'
QuestionCluster <- subset(Prostate_Low_Grade_Normal_Final, idents = "3")

###Merging data 
Lung_Complete_Clean <- merge(x=LungC_Clean_Final, y=LungN_Clean_Final)
View(Lung_Cancer_Final@meta.data)

TEST <- Lung_Complete_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

Lung_Complete_Clean <- TEST

DimPlot(Lung_Complete_Clean, reduction = "umap", label = TRUE, group.by="Condition")
VlnPlot(Lung_Complete_Clean, features = Tcellgenes)

#Check how many cells per patient per condition
LungN_Final$cellnumber <- paste(LungN_Final$PatientID, LungN_Final$Condition, sep = "_")
LungN_Final@meta.data$cellnumber %>% table()

saveRDS(Prostate_High_Grade_Complete_PatientsRemoved, file = "Hirz et al., 2023 (High-grade Prostate, n=3)")
#Genes used ---------- 
Tcellgenes <- c("CD3D", "CD8A", "CD3G", "CD3E", "CD8B", "LCK", "CD52")

NKcellgenes <- c("KLRF1", "GNLY", "NKG7", "PRF1")

NK1_module <- list(c('CD160','CTSD','CCL4','ADGRG1','CD38','CD247','CHST2','CX3CR1','KLRB1','LAIR2','IGFBP7','AKR1C3','FGFBP2','MYOM2','CLIC3','GZMB','PRF1','FCER1G','NKG7','SPON2'))

NK2_module <- list(c('LTB','FOS','IL2RB','IFITM3','COTL1','IL7R','PIK3R1','AREG','ZFP36L2','DUSP3','CD44','SELL','GPR183','CMC1','KLRC1','TCF7','TPT1','XCL2','XCL1','GZMK'))

NK3_module <- list(c('HBA1','CD3G','CD3D','DSTN','ZBTB38','CD2','PPDPF','PRDM1','LGALS1','S100A4','ITGB1','IL32','VIM','LINC01871','S100A5','PTMS','GZMH','CD3E','CCL5','KLRC2'))

ERsurvivalmarkers <- list(c("ERN1", "XBP1", "EIF2AK3", "ATF4", "ATF6", "HSPA5"))

ERstressmarkers <- list(c("HSPA1B", "CD52", "ATF3", "DUSP4", "JUN", "DNAJB1", "HSPA6", "RAB3GAP1", "HSPH1", "DNAJB4", "HSPB1", "HSPD1", "PPP1R14B", "GADD45G"))

Cytotoxicity <- list(c("CTSW", "PRF1", "GNLY", "GZMK", "GZMM", "GZMH", "GZMB", "GZMA"))

Inflammatory <- list(c("IL18", "IL15", "IL7", "IL6", "IL1B", "CXCL9", "CXCL10", "CCL5", "CCL4", "CCL3", "CCL2"))

#Ovarian Cancer ----------

###Reading in RDS file
Ovarian_normal_lymphoid_final <-  readRDS("~/Documents/areeba/Codes for analysis/Ovarian normal lymphoid final.rds")
Ovarian_tumor_lymphoid_final <-  readRDS("~/Documents/areeba/Codes for analysis/Ovarian tumor lymphoid final.rds")

###Changing assay to originalexp
Ovarian_normal_lymphoid_final[['originalexp']] = Ovarian_normal_lymphoid_final[['RNA']]
Ovarian_normal_lymphoid_final[['RNA']] = NULL
Ovarian_tumor_lymphoid_final[['originalexp']] = Ovarian_tumor_lymphoid_final[['RNA']]
Ovarian_tumor_lymphoid_final[['RNA']] = NULL

###Subsetting NK cells
#First check what cell types are present
Ovarian_normal_lymphoid_final@meta.data$celltype %>% table()
Ovarian_tumor_lymphoid_final@meta.data$celltype %>% table()
#Grab anything labeled NK except for NK-T and put into "idents"
Idents(Ovarian_normal_lymphoid_final) <- "celltype"
Idents(Ovarian_tumor_lymphoid_final) <- "celltype"
OvarianN = subset(Ovarian_normal_lymphoid_final, idents = "NK")
OvarianC = subset(Ovarian_tumor_lymphoid_final, idents=c("Circulating NK", "Tissue-resident NK"))
#Check that all cells have been subsetted
OvarianN@meta.data$celltype %>% table()
OvarianC@meta.data$celltype %>% table()

###Editing metadata 
#Edit metadata to have columns we want in the right order 
#(i.e. "nCount_originalexp", "nFeature_originalexp", "cancer_type", "celltype", "Patient", "Condition","Source")
##Ovarian Normal 
#Check the metadata to see what is present and what you need to add
View(OvarianN@meta.data)

#Editing to have common column names  
colnames(OvarianN@meta.data)
colnames(OvarianN@meta.data)[6] <- 'PatientID'
colnames(OvarianN@meta.data)[2] <- 'nCount_originalexp'
colnames(OvarianN@meta.data)[3] <- 'nFeature_originalexp'
colnames(OvarianN@meta.data)

#Changing patient names 
OvarianN@meta.data$PatientID %>% table()
Idents(OvarianN) <- 'PatientID'
OvarianN = RenameIdents(OvarianN, "Patient_11" = "Ovarian_001")
OvarianN$PatientID = Idents(OvarianN)
OvarianN = RenameIdents(OvarianN, "Patient_15" = "Ovarian_002")
OvarianN$PatientID = Idents(OvarianN)
Idents(OvarianN) <- 'PatientID'
OvarianN@meta.data$PatientID %>% table()
View(OvarianN)

#Adding missing columns to the metadata
OvarianN_Test <- addcolumns(OvarianN,"Ovarian","NK","Normal","Qian et al., 2020 (Ovarian)", "PatientID")
View(OvarianN_Test@meta.data)
OvarianN_Final <- OvarianN_Test
View(OvarianN_Final@meta.data)

#Running PCA 
TEST <- OvarianN_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:15)
TEST <- FindClusters(TEST, resolution = 1.2)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:15, reduction = "pca")

OvarianN_Final <- TEST

DimPlot(OvarianN_Final, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(OvarianN_Final, features = Tcellgenes)

Idents(OvarianN_Final) <- 'seurat_clusters'
cluster1 = subset(OvarianN_Final, idents = 1)
cluster0 = subset(cluster0, idents = 0)
OvarianN_Final = subset(OvarianN_Final, idents = c(0,2))
OvarianN_Final_Clean = subset(OvarianN_Final, idents = c(1,2,3,4,5))

##Ovarian Cancer
#Check the metadata to see what is present and what you need to add
View(OvarianC@meta.data)

#Editing to have common column names 
colnames(OvarianC@meta.data)
colnames(OvarianC@meta.data)[4] <- 'PatientID'
colnames(OvarianC@meta.data)[2] <- 'nCount_originalexp'
colnames(OvarianC@meta.data)[3] <- 'nFeature_originalexp'
colnames(OvarianC@meta.data)

#Changing patient names 
OvarianC@meta.data$PatientID %>% table()
Idents(OvarianC) <- 'PatientID'
OvarianC = RenameIdents(OvarianC, "Patient_11" = "Ovarian_001")
OvarianC$PatientID = Idents(OvarianC)
OvarianC = RenameIdents(OvarianC, "Patient_15" = "Ovarian_002")
OvarianC$PatientID = Idents(OvarianC)
Idents(OvarianC) <- 'PatientID'
OvarianC@meta.data$PatientID %>% table()

#Adding missing columns to the metadata
OvarianC_Test <- addcolumns(OvarianC,"Ovarian","NK","Tumor","Qian et al., 2020 (Ovarian)", "PatientID")
View(OvarianC_Test@meta.data)
OvarianC_Final <- OvarianC_Test
View(OvarianC_Final@meta.data)

#Running PCA
TEST <- OvarianC_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:16)
TEST <- FindClusters(TEST, resolution = 1.0)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:16, reduction = "pca")

OvarianC_Final <- TEST

DimPlot(OvarianC_Final, reduction = "umap", label = TRUE)

#Checking if there are a large significant population of T cells and filtering the data to remove these large populations 
VlnPlot(OvarianC_Final, features = Tcellgenes)

###Merging normal and tumor
Ovarian_Complete_Clean <- merge(x=OvarianN_Final, y=OvarianC_Final)
View(Ovarian_Complete_Clean@meta.data)

#Running PCA 
TEST <- Ovarian_Complete_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:11)
TEST <- FindClusters(TEST, resolution = 1.0)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:11, reduction = "pca")

Ovarian_Complete_Clean <- TEST

DimPlot(Ovarian_Complete_Clean, reduction = "umap", label = TRUE, group.by="Condition")

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(Ovarian_Complete_Clean, features = Tcellgenes)

#Check how many cells per patient per condition
Ovarian_Complete_Clean$cellnumber <- paste(Ovarian_Complete_Clean$PatientID, Ovarian_Complete_Clean$Condition, sep = "_")
Ovarian_Complete_Clean@meta.data$cellnumber %>% table()

saveRDS(Ovarian_Complete_Clean, file = "Qian et al., 2020 (Ovarian, n=2).rds")

#Colon Cancer ---------- 

###Reading in RDS file
MMRD_Final <- readRDS("~/Documents/areeba/Codes for analysis/MMRD_Tumor_Updated.rds")
MMRP_Final <- readRDS("~/Documents/areeba/Codes for analysis/MMRP_Tumor_Updated.rds")
Colon_Normal_Final <- readRDS("~/Documents/areeba/Codes for analysis/Colon_Normal_Updated.rds")
Colon_Normal_Final@meta.data$celltype %>% table()

###Subsetting NK cells
#First check what cell types are present
MMRP_Final@meta.data$celltype %>% table()
MMRD_Final@meta.data$celltype %>% table()
##Grab anything labeled NK except for NK-T and put into "idents"
Idents(Colon_Normal_Final) <- "celltype"
Idents(MMRP_Final) <- "celltype"
Idents(MMRD_Final) <- "celltype"
MMRD = subset(MMRD_Final, idents='NK')
MMRP = subset(MMRP_Final, idents= c('trNK', 'circulating NK'))
ColonN = subset(Colon_Normal_Final, idents='NK')
#Check that all cells have been subsetted
MMRP@meta.data$celltype %>% table()
MMRD@meta.data$celltype %>% table()
ColonN@meta.data$celltype %>% table()

###Editing the metadata 
#Edit metadata to have columns we want in the right order 
#(i.e. "nCount_originalexp", "nFeature_originalexp", "cancer_type", "celltype", "Patient", "Condition","Source")
##MMRD 
#Check the metadata to see what is present and what you need to add
View(MMRD@meta.data)

#Editing to have common column names 
colnames(MMRD@meta.data)
colnames(MMRD@meta.data)[8] <- 'PatientID'

#Changing patient names 
MMRD@meta.data$PatientID %>% table()
Idents(MMRD) <- 'PatientID'
MMRD = RenameIdents(MMRD, "Patient_C173_T" = "MMRD_001")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C170_T" = "MMRD_002")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C169_T" = "MMRD_003")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C168_T" = "MMRD_004")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C167_T" = "MMRD_005")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C165_T" = "MMRD_006")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C164_T" = "MMRD_007")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C163_T" = "MMRD_008")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C158_T" = "MMRD_009")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C156_T" = "MMRD_010")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C154_T" = "MMRD_011")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C152_T" = "MMRD_012")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C151_T" = "MMRD_013")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C147_T" = "MMRD_012")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C146_T" = "MMRD_013")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C144_T" = "MMRD_014")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C143_T" = "MMRD_015")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C142_T" = "MMRD_016")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C139_T" = "MMRD_017")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C138_T" = "MMRD_018")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C137_T" = "MMRD_019")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C132_T" = "MMRD_020")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C130_TA" = "MMRD_021A")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C130_TB" = "MMRD_021B")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C123_T" = "MMRD_022")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C122_T" = "MMRD_023")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C119_T" = "MMRD_024")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C118_T" = "MMRD_025")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C116_T" = "MMRD_026")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C115_T" = "MMRD_027")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C114_T" = "MMRD_028")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C111_T" = "MMRD_029")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C110_T" = "MMRD_030")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C109_T" = "MMRD_031")
MMRD$PatientID = Idents(MMRD)
MMRD = RenameIdents(MMRD, "Patient_C106_T" = "MMRD_032")
MMRD$PatientID = Idents(MMRD)
MMRD@meta.data$PatientID %>% table()

#Subsetting patients that have matching normal and tumor samples
Idents(MMRD) <- "PatientID"
MMRD = subset(MMRD, idents=c('MMRD_002', 'MMRD_013', 'MMRD_015', 'MMRD_016', 'MMRD_017', 'MMRD_018', 'MMRD_019', 'MMRD_020', 'MMRD_021A', 'MMRD_021B','MMRD_022', 'MMRD_023', 'MMRD_026', 'MMRD_027', 'MMRD_028', 'MMRD_029', 'MMRD_031', 'MMRD_032'))

#Adding missing columns to the metadata
MMRD_Test <- addcolumns(MMRD,"Colon","NK","Tumor","Pelka et al., 2021 (MMRD and MMRP)")
View(MMRD_Test@meta.data)
MMRD_Final <- MMRD_Test
View(MMRD_Final@meta.data)

#Runnning PCA
TEST <- MMRD_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:20)
TEST <- FindClusters(TEST, resolution = 1.0)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:20, reduction = "pca")

MMRD_Clean <- TEST

DimPlot(MMRD_Clean, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(MMRD_Clean, features = Tcellgenes)

Idents(MMRD_Final) <- 'seurat_clusters'
cluster7 = subset(MMRD_Final, idents = 7)
cluster7_clean = subset(cluster7_clean, idents = c(0,1,3,5))
cluster1 = subset(MMRD_Final, idents = 1)
cluster1_clean = subset(cluster1_clean, idents = c(0,1,3,5))
cluster10 = subset(MMRD_Final, idents = 10)
cluster10_clean = subset(cluster10, idents = 0)

MMRD_Clean <- subset(MMRD_Final, idents = c(0,2,3,4,5,6,8,9,11))
MMRD_Clean <- merge(x = MMRD_Clean, y = c(cluster1_clean, cluster10_clean, cluster7_clean))

##MMRP
#Check the metadata to see what is present and what you need to add
View(MMRP@meta.data)

#Editing to have common column names 
colnames(MMRP@meta.data)
colnames(MMRP@meta.data)[8] <- 'PatientID'

#Changing patient names 
MMRP@meta.data$PatientID %>% table()
Idents(MMRP) <- 'PatientID'
MMRP = RenameIdents(MMRP, "Patient_C172_T" = "MMRP_001")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C171_TA" = "MMRP_002")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C166_T" = "MMRP_003")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C162_T" = "MMRP_004")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C161_T" = "MMRP_005")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C160_T" = "MMRP_006")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C159_T" = "MMRP_007")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C157_T" = "MMRP_008")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C155_T" = "MMRP_009")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C153_T" = "MMRP_010")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C150_T" = "MMRP_011")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C149_T" = "MMRP_012")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C145_T" = "MMRP_013")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C140_T" = "MMRP_014")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C136_T" = "MMRP_015")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C135_T" = "MMRP_016")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C134_T" = "MMRP_017")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C133_T" = "MMRP_018")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C126_T" = "MMRP_019")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C125_T" = "MMRP_020")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C124_T" = "MMRP_021")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C113_T" = "MMRP_022")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C112_T" = "MMRP_023")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C107_T" = "MMRP_024")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C105" = "MMRP_025")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C104" = "MMRP_026")
MMRP$PatientID = Idents(MMRP)
MMRP = RenameIdents(MMRP, "Patient_C103" = "MMRP_027")
MMRP$PatientID = Idents(MMRP)
MMRP@meta.data$PatientID %>% table()

#Subsetting patients that have matching normal and tumor samples
Idents(MMRP) <- "PatientID"
MMRP = subset(MMRP, idents=c('MMRP_004','MMRP_008', 'MMRP_009', 'MMRP_014', 'MMRP_015', 'MMRP_016', 'MMRP_017', 'MMRP_018', 'MMRP_019', 'MMRP_020', 'MMRP_021', 'MMRP_022', 'MMRP_023', 'MMRP_024'))

#Adding missing columns to the metadata
MMRP_Test <- addcolumns(MMRP,"Colon","NK","Tumor","Pelka et al., 2021 (MMRD and MMRP)")
View(MMRP_Test@meta.data)
MMRP_Final <- MMRP_Test
View(MMRP_Final@meta.data)

#Running PCA
TEST <- MMRP_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 1.0)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

MMRP_Final <- TEST

DimPlot(MMRP_Final, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(MMRP_Final, features = Tcellgenes)

Idents(MMRP_Final) <- 'seurat_clusters'
QuestionCluster <- subset(MMRP_Final, idents = 2)
cluster4=FindMarkers(QuestionCluster, ident.1 = 4, min.pct = 0.25)

MMRP_Clean <- subset(MMRP_Final, idents = c(0,1,3))

##Colon Normal 
#Check the metadata to see what is present and what you need to add
View(ColonN@meta.data)

#Editing to have common column names 
colnames(ColonN@meta.data)
colnames(ColonN@meta.data)[8] <- 'PatientID'
ColonN@meta.data$PatientID %>% table()

#Changing patients name
Idents(ColonN) <- 'PatientID'
ColonN = RenameIdents(ColonN, "Patient_C170_N" = "Colon_001")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C165_N" = "MMRD_002")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C162_N" = "MMRP_004")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C157_N" = "MMRP_008")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C155_N" = "MMRP_009")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C151_N" = "MMRD_013")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C143_N" = "MMRD_015")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C142_N" = "MMRD_016")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C140_N" = "MMRP_014")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C139_N" = "MMRD_017")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C138_N" = "MMRD_018")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C137_N" = "MMRD_019")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C136_N" = "MMRP_015")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C135_N" = "MMRP_016")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C134_N" = "MMRP_017")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C133_N" = "MMRP_018")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C132_N" = "MMRD_020")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C130_N" = "MMRD_021")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C129_N" = "Colon_002")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C126_N" = "MMRP_019")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C125_N" = "MMRP_020")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C124_N" = "MMRP_021")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C123_N" = "MMRD_022")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C122_N" = "MMRD_023")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C116_N" = "MMRD_026")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C115_N" = "MMRD_027")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C114_N" = "MMRD_028")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C113_N" = "MMRP_022")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C112_N" = "MMRP_023")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C111_N" = "MMRD_029")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C110_N" = "MMRD_030")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C109_N" = "MMRD_031")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C107_N" = "MMRP_024")
ColonN$PatientID = Idents(ColonN)
ColonN = RenameIdents(ColonN, "Patient_C106_N" = "MMRD_032")
ColonN$PatientID = Idents(ColonN)
ColonN@meta.data$PatientID %>% table()

#Subsetting patients that have matching normal and tumor samples
Idents(ColonN) <- "PatientID"
ColonN = subset(ColonN, idents=c('MMRP_004','MMRP_008', 'MMRP_009', 'MMRP_014', 'MMRP_015', 'MMRP_016', 'MMRP_017', 'MMRP_018', 'MMRP_019', 'MMRP_020', 'MMRP_021', 'MMRP_022', 'MMRP_023', 'MMRP_024', 'MMRD_002', 'MMRD_013', 'MMRD_015', 'MMRD_016', 'MMRD_017', 'MMRD_018', 'MMRD_019', 'MMRD_020', 'MMRD_021','MMRD_022', 'MMRD_023', 'MMRD_026', 'MMRD_027', 'MMRD_028', 'MMRD_029', 'MMRD_031', 'MMRD_032'))

#Adding missing columns to the metadata
ColonN_Test <- addcolumns(ColonN,"Colon","NK","Normal","Pelka et al., 2021 (MMRD and MMRP)")
View(ColonN_Test@meta.data)
ColonN_Final <- ColonN_Test
View(ColonN_Final@meta.data)

#Running PCA
TEST <- ColonN_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

ColonN_Final <- TEST

DimPlot(ColonN_Final, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(ColonN_Final, features = Tcellgenes)

ColonN_Clean <- ColonN_Final


###Merging normal and tumor 
Colon_Complete_Clean <- merge(x=ColonN_Clean, y=c(MMRD_Clean,MMRP_Clean))
View(Colon_Complete_Clean@meta.data)

#Running PCA
TEST <- Colon_Complete_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:11)
TEST <- FindClusters(TEST, resolution = 1.0)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:11, reduction = "pca")

Colon_Complete_Clean <- TEST

DimPlot(Colon_Complete_Clean, reduction = "umap", label = TRUE, group.by="Condition")

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(Colon_Complete_Clean, features = Tcellgenes)

cluster12 = subset(Colon_Complete_Clean, idents = 12)
cluster2_clean = subset(cluster2, idents =c(0,1,2,3))

Colon_Complete_Clean = subset(Colon_Complete_Clean, idents = c(0,1,2,3,4,5,6,8,9,10,12,13))
Colon_Complete_Clean = merge(x = Colon_Complete_Clean, y = cluster2_clean)

#Check how many cells per patient per condition
Colon_Complete_Clean$cellnumber <- paste(Colon_Complete_Clean$PatientID, Colon_Complete_Clean$Condition, sep = "_")
Colon_Complete_Clean@meta.data$Condition %>% table()

saveRDS(Colon_Complete_PatientsRemoved, file = "Pelka et al., 2021 (MMRD, n=3 and MMRP, n=2)")
#Kidney Cancer ---------- 

###Reading in RDS file
Kidney_Cancer_Final <- readRDS("~/Documents/areeba/Codes for analysis/Kidney tumor final.rds")
Kidney_Normal_Final <- readRDS("~/Documents/areeba/Codes for analysis/Kidney normal final.rds")

###Changing assay to originalexp
Kidney_Cancer_Final[['originalexp']] = Kidney_Cancer_Final[['RNA']]
Kidney_Cancer_Final[['RNA']] = NULL
Kidney_Normal_Final[['originalexp']] = Kidney_Normal_Final[['RNA']]
Kidney_Normal_Final[['RNA']] = NULL

###Subsetting NK cells
#First check what cell types are present
Kidney_Cancer_Final@meta.data$celltype %>% table()
Kidney_Normal_Final@meta.data$celltype %>% table()
##Grab anything labeled NK except for NK-T and put into "idents"
Idents(Kidney_Cancer_Final) <- "celltype"
Idents(Kidney_Normal_Final) <- "celltype"
KidneyC = subset(Kidney_Cancer_Final, idents=c('Tissue-resident NK', 'Circulating NK'))
View(KidneyC@meta.data)
KidneyN = subset(Kidney_Normal_Final, idents=c("NK", "Tissue-resident NK"))
View(KidneyN@meta.data)
#Check that all of the cells have been subsetted
KidneyC@meta.data$celltype %>% table()
KidneyN@meta.data$celltype %>% table()

###Editing the metadata
#Edit metadata to have columns we want in the right order 
#(i.e. "nCount_originalexp", "nFeature_originalexp", "cancer_type", "celltype", "Patient", "Condition","Source")
##Kidney Cancer
#Check the metadata to see what is present and what you need to add
View(KidneyC@meta.data)

#Editing to have common column names 
colnames(KidneyC@meta.data)
colnames(KidneyC@meta.data)[18] <- 'PatientID'
colnames(KidneyC@meta.data)[2] <- 'nCount_originalexp'
colnames(KidneyC@meta.data)[3] <- 'nFeature_originalexp'
colnames(KidneyC@meta.data)

#Changing patient names 
KidneyC@meta.data$PatientID %>% table()
Idents(KidneyC) <- 'PatientID'
KidneyC = RenameIdents(KidneyC, "p022" = "Kidney_001")
KidneyC$PatientID = Idents(KidneyC)
KidneyC = RenameIdents(KidneyC, "p027" = "Kidney_002")
KidneyC$PatientID = Idents(KidneyC)
KidneyC = RenameIdents(KidneyC, "p029" = "Kidney_003")
KidneyC$PatientID = Idents(KidneyC)
KidneyC@meta.data$PatientID %>% table()

#Adding missing columns to the metadata
KidneyC_Test <- addcolumns(KidneyC,"Kidney","NK","Tumor","Massenet-Regad et al., 2023 (ccRCC)", "PatientID")
View(KidneyC_Test@meta.data)
KidneyC_Final <- KidneyC_Test
View(KidneyC_Final@meta.data)

#Running PCA
TEST <- KidneyC_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 1.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

KidneyC_Final <- TEST

DimPlot(QuestionCluster, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(QuestionCluster, features = Tcellgenes)

Idents(KidneyC_Final) <- 'seurat_clusters'
QuestionCluster <- subset(KidneyC_Final, idents = 5)
cluster4_clean <- subset(QuestionCluster, idents = c(0,3,4,5))
cluster5=FindMarkers(QuestionCluster, ident.1 = 5, min.pct = 0.25)
KidneyC_Final_Clean <- subset(KidneyC_Final, idents = c(0,1,2,3,4,6,7))
KidneyC_Clean <- merge(x = cluster4_clean, y = KidneyC_Final_Clean)
View(KidneyC_Clean@meta.data)

##Kidney Normal 
#Check the metadata to see what is present and what you need to add
View(KidneyN@meta.data)

#Editing to have common column names 
colnames(KidneyN@meta.data)
colnames(KidneyN@meta.data)[18] <- 'PatientID'
colnames(KidneyN@meta.data)[2] <- 'nCount_originalexp'
colnames(KidneyN@meta.data)[3] <- 'nFeature_originalexp'

#Changing patient names 
KidneyN@meta.data$PatientID %>% table()
Idents(KidneyN) <- 'PatientID'
KidneyN = RenameIdents(KidneyN, "p022" = "Kidney_001")
KidneyN$PatientID = Idents(KidneyN)
KidneyN = RenameIdents(KidneyN, "p027" = "Kidney_002")
KidneyN$PatientID = Idents(KidneyN)
KidneyN = RenameIdents(KidneyN, "p029" = "Kidney_003")
KidneyN$PatientID = Idents(KidneyN)
KidneyN@meta.data$PatientID %>% table()

#Adding missing columns to the metadata
KidneyN_Test <- addcolumns(KidneyN,"Kidney","NK","Normal","Massenet-Regad et al., 2023 (ccRCC)")
View(KidneyN_Test@meta.data)
KidneyN_Final <- KidneyN_Test
View(KidneyN_Final@meta.data)

#Running PCA
TEST <- KidneyN_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

KidneyN_Final <- TEST

DimPlot(KidneyN_Final, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(KidneyN_Final, features = Tcellgenes)

KidneyN_Clean <- KidneyN_Final
View(KidneyN_Clean@meta.data)

###Merging normal and tumor
Kidney_Complete_Clean <- merge(x=KidneyC_Clean, y=KidneyN_Clean)
View(Kidney_Complete_Clean@meta.data)

#Running PCA
TEST <- Kidney_Complete_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:11)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:11, reduction = "pca")

Kidney_Complete_Clean <- TEST

DimPlot(Kidney_Complete_Clean, reduction = "umap", label = TRUE, group.by="Condition")

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(Kidney_Complete_Clean, features = Tcellgenes)

#Check how many cells per patient per condition
Kidney_Complete_Clean$cellnumber <- paste(Kidney_Complete_Clean$PatientID, Kidney_Complete_Clean$Condition, sep = "_")
Kidney_Complete_Clean@meta.data$cellnumber %>% table()

saveRDS(Kidney_Complete_Clean, file = "Massenet-Regad et al., 2023 (ccRCC, n=3)")
#Lung Cancer ----------
###Reading in RDS file
Lung_Cancer_Final <- readRDS("~/Documents/areeba/Codes for analysis/Lung_Tumor.rds")
Lung_Normal_Final <- readRDS("~/Documents/areeba/Codes for analysis/Lung_Normal.rds")

###Subsetting NK cells
#First check what cell types are present
Lung_Cancer_Final@meta.data$celltype %>% table()
Lung_Normal_Final@meta.data$celltype %>% table()
Idents(Lung_Cancer_Final) <- "celltype"
Idents(Lung_Normal_Final) <- "celltype"
#Grab anything labeled NK except for NK-T and put into "idents"
LungC = subset(Lung_Cancer_Final, idents='NK')
LungN = subset(Lung_Normal_Final, idents='NK')
#Check that all cells have been subsetted
LungC@meta.data$celltype %>% table()
LungN@meta.data$celltype %>% table()

###Editing metadata
#Edit metadata to have columns we want in the right order 
#(i.e. "nCount_originalexp", "nFeature_originalexp", "cancer_type", "celltype", "Patient", "Condition","Source")
##Lung Normal
#Check the metadata to see what is present and what you need to add
View(LungN@meta.data)

#Editing to have common column names 
colnames(LungN@meta.data)
colnames(LungN@meta.data)[9] <- 'PatientID'

#Changing patient names
LungN@meta.data$PatientID %>% table()
Idents(LungN) <- 'PatientID'
LungN$PatientID = Idents(LungN)
LungN = RenameIdents(LungN, "Patient_BT1301" = "Lung_005")
LungN$PatientID = Idents(LungN)
LungN = RenameIdents(LungN, "Patient_BT1294" = "Lung_004")
LungN$PatientID = Idents(LungN)
LungN = RenameIdents(LungN, "Patient_BT1293" = "Lung_003")
LungN$PatientID = Idents(LungN)
LungN = RenameIdents(LungN, "Patient_BT1247" = "Lung_002")
LungN$PatientID = Idents(LungN)
LungN@meta.data$PatientID %>% table()

#Adding missing columns to the metadata
LungN_Test <- addcolumns(LungN,"Lung","NK","Normal","Lambrechts et al., 2018 (Lung)", "PatientID")
View(LungN_Final@meta.data)
LungN_Final <- LungN_Test
View(LungN_Final@meta.data)

#Running PCA
TEST <- LungN_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:11)
TEST <- FindClusters(TEST, resolution = 1.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

LungN_Final <- TEST

DimPlot(LungN_Final, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(LungN_Final, features = NKcellgenes)

Idents(LungN_Final) <- 'seurat_clusters'
QuestionCluster <- subset(LungN_Final, idents = 0)
cluster0_clean <- subset(QuestionCluster, idents=c(0,1,2,3,6,7,8,9))
LungN_Clean_Final <- subset(LungN_Final, idents = c(1,2,3))
LungN_Clean <- merge(x = cluster0_clean, y = LungC_Clean_Final)

##Lung Cancer 
#Check the metadata to see what is present and what you need to add
View(LungC@meta.data)

#Editing to have common column names 
colnames(LungC@meta.data)
colnames(LungC@meta.data)[9] <- 'PatientID'

#Changing patient names 
LungC@meta.data$PatientID %>% table()
Idents(LungC) <- 'PatientID'
LungC = RenameIdents(LungC, "Patient_BT2B" = "Lung_002")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT2A" = "Lung_002")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT1C" = "Lung_006")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT1B" = "Lung_007")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT1A" = "Lung_008")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT1300" = "Lung_005")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT1299" = "Lung_005")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT1298" = "Lung_005")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT1297" = "Lung_004")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT1296" = "Lung_004")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT1295" = "Lung_004")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT1292" = "Lung_003")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT1291" = "Lung_003")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT1290" = "Lung_003")
LungC$PatientID = Idents(LungC)
LungC = RenameIdents(LungC, "Patient_BT1249" = "Lung_002")
LungC$PatientID = Idents(LungC)

#Subsetting patients that have matching normal and tumor samples
LungC@meta.data$PatientID %>% table()
LungC = subset(LungC, idents=c('Lung_002', 'Lung_003', 'Lung_004', 'Lung_005'))
View(LungC)

#Adding missing columns to the metadata
LungC_Test <- addcolumns(LungC,"Lung","NK","Tumor","Lambrechts et al., 2018 (Lung)", "PatientID")
View(LungC_Test@meta.data)
LungC_Final <- LungC_Test
View(LungC_Final@meta.data)

#Running PCA
TEST <- LungC_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:11)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:11, reduction = "pca")

LungC_Final <- TEST

DimPlot(LungC_Final, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(LungC_Final, features = Tcellgenes)

Idents(LungC_Final) <- 'seurat_clusters'
QuestionCluster <- subset(LungC_Final, idents = 0)
cluster0 <- FindMarkers(LungC_Final, ident.1=0)

LungC_Clean <- subset(LungC_Final, idents = c(1,2,3))

###Merging normal and tumor 
Lung_Complete_Clean <- merge(x=LungC_Clean, y=LungN_Clean)
View(Lung_Complete_Clean@meta.data)

#Running PCA
TEST <- Lung_Complete_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

Lung_Complete_Clean <- TEST

DimPlot(Lung_Complete_Clean, reduction = "umap", label = TRUE, group.by="Condition")

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(Lung_Complete_Clean, features = Tcellgenes)

#Check how many cells per patient per condition
Lung_Complete_Clean$cellnumber <- paste(Lung_Complete_Clean$PatientID, Lung_Complete_Clean$Condition, sep = "_")
Lung_Complete_Clean@meta.data$Condition %>% table()

saveRDS(Lung_Complete_PatientsRemoved, file = "Lambrechts et al., 2018 (Lung, n=1)")
#Brain Cancer ----------

###Reading in RDS file
Brain_Combined <- readRDS("~/Documents/areeba/Codes for analysis/Combined_Glioblastoma.rds")
Brain_Combined@meta.data$celltype %>% table()

###Changing assay to originalexp
Brain_Combined[['originalexp']] = Brain_Combined[['RNA']]
Brain_Combined[['RNA']] = NULL

###Subsetting NK cells
Brain_Combined@meta.data$celltype %>% table()
Idents(Brain_Combined) <- "celltype"
#Grab anything labeled NK except for NK-T and put into "idents"
Brain_Combined = subset(Brain_Combined, idents= 'NK cells')
#Check that all cells have been subsetted
Brain_Combined@meta.data$celltype %>% table()

###Editing the metadata
#Edit metadata to have columns we want in the right order 
#(i.e. "nCount_originalexp", "nFeature_originalexp", "cancer_type", "celltype", "Patient", "Condition","Source")
##Brain Cancer
#Check the metadata to see what is present and what you need to add
View(Brain_Combined@meta.data)

#Editing to have common column names 
colnames(Brain_Combined@meta.data)
colnames(Brain_Combined@meta.data)[5] <- 'Condition'
colnames(Brain_Combined@meta.data)[17] <- 'PatientID'
colnames(Brain_Combined@meta.data)[2] <- 'nCount_originalexp'
colnames(Brain_Combined@meta.data)[3] <- 'nFeature_originalexp'
colnames(Brain_Combined@meta.data)

#Subsetting Cancer NK cells
Idents(Brain_Combined) <- "Condition"
BrainC = subset(Brain_Combined, idents = 'glioblastoma')
View(BrainC@meta.data)
BrainC@meta.data$PatientID %>% table()

#Changing patient names 
Idents(BrainC) <- 'PatientID'
BrainC$PatientID = Idents(BrainC)
BrainC = RenameIdents(BrainC, "patient 604" = "Brain_001")
BrainC$PatientID = Idents(BrainC)
BrainC = RenameIdents(BrainC, "patient 622" = "Brain_002")
BrainC$PatientID = Idents(BrainC)
BrainC = RenameIdents(BrainC, "patient 617" = "Brain_003")
BrainC$PatientID = Idents(BrainC)
BrainC = RenameIdents(BrainC, "patient 616" = "Brain_004")
BrainC$PatientID = Idents(BrainC)
View(BrainC@meta.data)

#Adding missing columns to the metadata
BrainC_Test <- addcolumns(BrainC,"Brain","NK","Tumor","Schmassmann et al., 2023 (GBM)", "PatientID")
View(BrainC_Final@meta.data)
BrainC_Final <- BrainC_Test

#Running PCA
TEST <- BrainC_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:11)
TEST <- FindClusters(TEST, resolution = 1.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

BrainC_Clean <- TEST

DimPlot(BrainC_Clean, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(BrainC_Clean, features = Tcellgenes)

Idents(BrainC_Final) <- 'seurat_clusters'
QuestionCluster <- subset(BrainC_Final, idents = 0)
BrainC_Final_Clean <- subset(BrainC_Final, idents = c(1,2,3))
BrainC_Clean <- BrainC_Final_Clean

##Brain Normal
#Check the metadata to see what is present and what you need to add
View(Brain_Combined@meta.data)

#Subsetting Normal NK cells 
Idents(Brain_Combined) <- "Condition"
BrainN = subset(Brain_Combined, idents = 'normal')
View(BrainN@meta.data)

#Changing patient names
BrainN@meta.data$PatientID %>% table()
Idents(BrainN) <- 'PatientID'
BrainN$PatientID = Idents(BrainN)
BrainN = RenameIdents(BrainN, "patient 604" = "Brain_001")
BrainN$PatientID = Idents(BrainN)
BrainN = RenameIdents(BrainN, "patient 622" = "Brain_002")
BrainN$PatientID = Idents(BrainN)
BrainN = RenameIdents(BrainN, "patient 617" = "Brain_003")
BrainN$PatientID = Idents(BrainN)
BrainN = RenameIdents(BrainN, "patient 616" = "Brain_004")
BrainN$PatientID = Idents(BrainN)
BrainN@meta.data$PatientID %>% table()
View(BrainN@meta.data)

#Adding missing columns to the metadata
BrainN_Test <- addcolumns(BrainN,"Brain","NK","Normal","Schmassmann et al., 2023 (GBM)", "PatientID")
View(BrainN_Test@meta.data)
BrainN_Final <- BrainN_Test
View(BrainN_Final@meta.data)

#Running PCA
TEST <- BrainN_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:11)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

BrainN_Final <- TEST

DimPlot(BrainN_Final, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations 
VlnPlot(BrainN_Final, features = Tcellgenes)

BrainN_Clean <- BrainN_Final

###Merging normal and tumor 
Brain_Complete_Clean <- merge(x=BrainN_Clean, y=BrainC_Clean)
View(Brain_Complete_Clean@meta.data)

#Running PCA
TEST <- Brain_Complete_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:11)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

Brain_Complete_Clean <- TEST

DimPlot(Brain_Complete_Clean, reduction = "umap", label = TRUE, group.by="Condition")

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(Brain_Complete_Clean, features = Tcellgenes)

#Check how many cells per patient per condition
Brain_Complete_Clean$cellnumber <- paste(Brain_Complete_Clean$PatientID, Brain_Complete_Clean$Condition, sep = "_")
Brain_Complete_Clean@meta.data$Condition %>% table()

saveRDS(Brain_Complete_Clean, file = "Schmassmann et al., 2023 (GBM, n=4)")
#Prostate Cancer ---------- 

###Reading in RDS file
Prostate_High_Grade_Normal <- readRDS("~/Documents/areeba/Codes for analysis/Prostate High Grade Normal.rds")
Prostate_High_Grade_Tumor <- readRDS("~/Documents/areeba/Codes for analysis/Prostate High Grade Tumor.rds")
Prostate_Low_Grade_Normal <- readRDS("~/Documents/areeba/Codes for analysis/Prostate Low Grade Normal.rds")
Prostate_Low_Grade_Tumor <- readRDS("~/Documents/areeba/Codes for analysis/Prostate Low Grade Tumor.rds")

###Changing assay to originalexp
Prostate_High_Grade_Normal[['originalexp']] = Prostate_High_Grade_Normal[['RNA']]
Prostate_High_Grade_Normal[['RNA']] = NULL
Prostate_High_Grade_Tumor[['originalexp']] = Prostate_High_Grade_Tumor[['RNA']]
Prostate_High_Grade_Tumor[['RNA']] = NULL
Prostate_Low_Grade_Normal[['originalexp']] = Prostate_Low_Grade_Normal[['RNA']]
Prostate_Low_Grade_Normal[['RNA']] = NULL
Prostate_Low_Grade_Tumor[['originalexp']] = Prostate_Low_Grade_Tumor[['RNA']]
Prostate_Low_Grade_Tumor[['RNA']] = NULL


###Subsetting NK cells
#First check what cell types are present
Prostate_High_Grade_Normal@meta.data$celltype %>% table()
Prostate_High_Grade_Normal@meta.data$celltype %>% table()
Prostate_Low_Grade_Normal@meta.data$celltype %>% table()
Prostate_Low_Grade_Tumor@meta.data$celltype %>% table()
#Grab anything labeled NK except for NK-T and put into "idents"
Idents(Prostate_Low_Grade_Tumor) <- "celltype"
Idents(Prostate_High_Grade_Tumor) <- "celltype"
Idents(Prostate_Low_Grade_Normal) <- "celltype"
Idents(Prostate_High_Grade_Normal) <- "celltype"
Prostate_High_Grade_Normal = subset(Prostate_High_Grade_Normal, idents=c('NK', 'Tissue-resident NK'))
Prostate_High_Grade_Tumor = subset(Prostate_High_Grade_Tumor, idents= 'NK')
Prostate_Low_Grade_Normal = subset(Prostate_Low_Grade_Normal, idents= 'NK')
Prostate_Low_Grade_Tumor = subset(Prostate_Low_Grade_Tumor, idents=c('Circulating NK', 'Tissue-resident NK'))
#Check that all cells have been subsetted
Prostate_High_Grade_Normal@meta.data$celltype %>% table()
Prostate_High_Grade_Tumor@meta.data$celltype %>% table()
Prostate_Low_Grade_Normal@meta.data$celltype %>% table()
Prostate_Low_Grade_Tumor@meta.data$celltype %>% table()

###Editing metadata
#Edit metadata to have columns we want in the right order 
#(i.e. "nCount_originalexp", "nFeature_originalexp", "cancer_type", "celltype", "Patient", "Condition","Source")
##Prostate High Grade Normal 
#Check the metadata to see what is present and what you need to add
View(Prostate_High_Grade_Normal)

#Editing to have common column names 
colnames(Prostate_High_Grade_Normal@meta.data)
colnames(Prostate_High_Grade_Normal@meta.data)[1] <- 'PatientID'
colnames(Prostate_High_Grade_Normal@meta.data)[2] <- 'nCount_originalexp'
colnames(Prostate_High_Grade_Normal@meta.data)[3] <- 'nFeature_originalexp'

#Changing patient names
Prostate_High_Grade_Normal@meta.data$PatientID %>% table()
Idents(Prostate_High_Grade_Normal) <- 'PatientID'
Prostate_High_Grade_Normal$PatientID = Idents(Prostate_High_Grade_Normal)
Prostate_High_Grade_Normal = RenameIdents(Prostate_High_Grade_Normal, "SCG.PCA15.N.HG" = "Prostate_High_Grade_001")
Prostate_High_Grade_Normal$PatientID = Idents(Prostate_High_Grade_Normal)
Prostate_High_Grade_Normal = RenameIdents(Prostate_High_Grade_Normal, "SCG.PCA19.N.HG" = "Prostate_High_Grade_002")
Prostate_High_Grade_Normal$PatientID = Idents(Prostate_High_Grade_Normal)
Prostate_High_Grade_Normal = RenameIdents(Prostate_High_Grade_Normal, "SCG.PCA22.N.HG" = "Prostate_High_Grade_003")
Prostate_High_Grade_Normal$PatientID = Idents(Prostate_High_Grade_Normal)
Prostate_High_Grade_Normal = RenameIdents(Prostate_High_Grade_Normal, "SCG.PCA4.N.HG" = "Prostate_High_Grade_004")
Prostate_High_Grade_Normal$PatientID = Idents(Prostate_High_Grade_Normal)
Prostate_High_Grade_Normal = RenameIdents(Prostate_High_Grade_Normal, "SCG.PCA6.N.HG" = "Prostate_High_Grade_005")
Prostate_High_Grade_Normal$PatientID = Idents(Prostate_High_Grade_Normal)
Prostate_High_Grade_Normal@meta.data$PatientID %>% table()

#Adding missing columns to the metadata
Prostate_High_Grade_Normal_Test <- addcolumns(Prostate_High_Grade_Normal,"Prostate","NK","Normal","Hirz et al., 2023 (Prostate)", "PatientID")
View(Prostate_High_Grade_Normal_Final@meta.data)
Prostate_High_Grade_Normal_Final <- Prostate_High_Grade_Normal_Test
View(Prostate_High_Grade_Normal_Final@meta.data)

#Running PCA
TEST <- Prostate_High_Grade_Normal_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

Prostate_High_Grade_Normal_Clean <- TEST

DimPlot(Prostate_High_Grade_Normal_Clean, reduction = "umap", label = TRUE, group.by = "Condition")

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(Prostate_High_Grade_Normal_Clean, features = Tcellgenes)

Prostate_High_Grade_Normal_Clean <- Prostate_High_Grade_Normal_Final

##Prostate High Grade Cancer 
##Check the metadata to see what is present and what you need to add
View(Prostate_High_Grade_Tumor)

#Editing to have common column names 
colnames(Prostate_High_Grade_Tumor@meta.data)
colnames(Prostate_High_Grade_Tumor@meta.data)[1] <- 'PatientID'
colnames(Prostate_High_Grade_Tumor@meta.data)[2] <- 'nCount_originalexp'
colnames(Prostate_High_Grade_Tumor@meta.data)[3] <- 'nFeature_originalexp'

#Changing patient names
Prostate_High_Grade_Tumor@meta.data$PatientID %>% table()
Idents(Prostate_High_Grade_Tumor) <- 'PatientID'
Prostate_High_Grade_Tumor$PatientID = Idents(Prostate_High_Grade_Tumor)
Prostate_High_Grade_Tumor = RenameIdents(Prostate_High_Grade_Tumor, "SCG.PCA15.T.HG" = "Prostate_High_Grade_001")
Prostate_High_Grade_Tumor$PatientID = Idents(Prostate_High_Grade_Tumor)
Prostate_High_Grade_Tumor = RenameIdents(Prostate_High_Grade_Tumor, "SCG.PCA19.T.HG" = "Prostate_High_Grade_002")
Prostate_High_Grade_Tumor$PatientID = Idents(Prostate_High_Grade_Tumor)
Prostate_High_Grade_Tumor = RenameIdents(Prostate_High_Grade_Tumor, "SCG.PCA22.T.HG" = "Prostate_High_Grade_003")
Prostate_High_Grade_Tumor$PatientID = Idents(Prostate_High_Grade_Tumor)
Prostate_High_Grade_Tumor = RenameIdents(Prostate_High_Grade_Tumor, "SCG.PCA4.T.HG" = "Prostate_High_Grade_004")
Prostate_High_Grade_Tumor$PatientID = Idents(Prostate_High_Grade_Tumor)
Prostate_High_Grade_Tumor = RenameIdents(Prostate_High_Grade_Tumor, "SCG.PCA6.T.HG" = "Prostate_High_Grade_005")
Prostate_High_Grade_Tumor$PatientID = Idents(Prostate_High_Grade_Tumor)
Prostate_High_Grade_Tumor@meta.data$PatientID %>% table()

#Adding missing columns to the metadata
Prostate_High_Grade_Tumor_Test <- addcolumns(Prostate_High_Grade_Tumor,"Prostate","NK","Tumor","Hirz et al., 2023 (Prostate)", "PatientID")
View(Prostate_High_Grade_Tumor_Final@meta.data)
Prostate_High_Grade_Tumor_Final <- Prostate_High_Grade_Tumor_Test
View(Prostate_High_Grade_Tumor_Final@meta.data)

#Running PCA
TEST <- Prostate_High_Grade_Tumor_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

Prostate_High_Grade_Tumor_Final <- TEST

DimPlot(Prostate_High_Grade_Tumor_Final, reduction = "umap", label = TRUE, group.by = "Condition")

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(Prostate_High_Grade_Tumor_Final, features = Tcellgenes)

Prostate_High_Grade_Tumor_Clean <- Prostate_High_Grade_Tumor_Final

###Merging normal and tumor 
Prostate_High_Grade_Complete_Clean <- merge(x=Prostate_High_Grade_Normal_Clean, y=Prostate_High_Grade_Tumor_Clean)
View(Prostate_High_Grade_Complete_Clean@meta.data)

#Running PCA
TEST <- Prostate_High_Grade_Complete_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

Prostate_High_Grade_Complete_Clean <- TEST

DimPlot(Prostate_High_Grade_Complete_Clean, reduction = "umap", label = TRUE, group.by="Condition")
VlnPlot(Prostate_High_Grade_Complete_Clean, features = Tcellgenes)
VlnPlot(Prostate_High_Grade_Complete_Clean, features = NKcellgenes)

#Check how many cells per patient per condition
Prostate_High_Grade_Complete_Clean$cellnumber <- paste(Prostate_High_Grade_Complete_Clean$PatientID, Prostate_High_Grade_Complete_Clean$Condition, sep = "_")
Prostate_High_Grade_Complete_Clean@meta.data$Condition %>% table()

saveRDS(Prostate_High_Grade_Complete_PatientsRemoved, file = "Hirz et al., 2023 (High-grade Prostate, n=3)")

##Prostate Low Grade Normal
#Check the metadata to see what is present and what you need to add
View(Prostate_Low_Grade_Normal)

#Editing to have common column names  
colnames(Prostate_Low_Grade_Normal@meta.data)
colnames(Prostate_Low_Grade_Normal@meta.data)[1] <- 'PatientID'
colnames(Prostate_Low_Grade_Normal@meta.data)[2] <- 'nCount_originalexp'
colnames(Prostate_Low_Grade_Normal@meta.data)[3] <- 'nFeature_originalexp'

#Changing patients names 
Prostate_Low_Grade_Normal@meta.data$PatientID %>% table()
Idents(Prostate_Low_Grade_Normal) <- 'PatientID'
Prostate_Low_Grade_Normal$PatientID = Idents(Prostate_Low_Grade_Normal)
Prostate_Low_Grade_Normal = RenameIdents(Prostate_Low_Grade_Normal, "SCG.PCA11.N.LG" = "Prostate_Low_Grade_001")
Prostate_Low_Grade_Normal$PatientID = Idents(Prostate_Low_Grade_Normal)
Prostate_Low_Grade_Normal = RenameIdents(Prostate_Low_Grade_Normal, "SCG.PCA12.N.LG" = "Prostate_Low_Grade_002")
Prostate_Low_Grade_Normal$PatientID = Idents(Prostate_Low_Grade_Normal)
Prostate_Low_Grade_Normal = RenameIdents(Prostate_Low_Grade_Normal, "SCG.PCA17.N.LG" = "Prostate_Low_Grade_003")
Prostate_Low_Grade_Normal$PatientID = Idents(Prostate_Low_Grade_Normal)
Prostate_Low_Grade_Normal = RenameIdents(Prostate_Low_Grade_Normal, "SCG.PCA18.N.LG" = "Prostate_Low_Grade_004")
Prostate_Low_Grade_Normal$PatientID = Idents(Prostate_Low_Grade_Normal)
Prostate_Low_Grade_Normal = RenameIdents(Prostate_Low_Grade_Normal, "SCG.PCA20.N.LG" = "Prostate_Low_Grade_005")
Prostate_Low_Grade_Normal$PatientID = Idents(Prostate_Low_Grade_Normal)
Prostate_Low_Grade_Normal = RenameIdents(Prostate_Low_Grade_Normal, "SCG.PCA21.N.LG" = "Prostate_Low_Grade_006")
Prostate_Low_Grade_Normal$PatientID = Idents(Prostate_Low_Grade_Normal)
Prostate_Low_Grade_Normal = RenameIdents(Prostate_Low_Grade_Normal, "SCG.PCA3.N.LG" = "Prostate_Low_Grade_007")
Prostate_Low_Grade_Normal$PatientID = Idents(Prostate_Low_Grade_Normal)
Prostate_Low_Grade_Normal = RenameIdents(Prostate_Low_Grade_Normal, "SCG.PCA5.N.LG" = "Prostate_Low_Grade_008")
Prostate_Low_Grade_Normal$PatientID = Idents(Prostate_Low_Grade_Normal)
Prostate_Low_Grade_Normal = RenameIdents(Prostate_Low_Grade_Normal, "SCG.PCA9.N.LG" = "Prostate_Low_Grade_009")
Prostate_Low_Grade_Normal$PatientID = Idents(Prostate_Low_Grade_Normal)
Prostate_Low_Grade_Normal@meta.data$PatientID %>% table()

#Subsetting patients that have matching normal and tumor samples
Idents(Prostate_Low_Grade_Normal) <- 'PatientID'
Prostate_Low_Grade_Normal = subset(Prostate_Low_Grade_Normal, idents=c('Prostate_Low_Grade_009', 'Prostate_Low_Grade_008', 'Prostate_Low_Grade_007', 'Prostate_Low_Grade_006', 'Prostate_Low_Grade_005', 'Prostate_Low_Grade_004', 'Prostate_Low_Grade_003', 'Prostate_Low_Grade_002'))
Prostate_Low_Grade_Normal@meta.data$PatientID %>% table()

#Adding missing columns to the metadata
Prostate_Low_Grade_Normal_Test <- addcolumns(Prostate_Low_Grade_Normal,"Prostate","NK","Normal","Hirz et al., 2023 (Prostate)", "PatientID")
View(Prostate_Low_Grade_Normal_Final@meta.data)
Prostate_Low_Grade_Normal_Final <- Prostate_Low_Grade_Normal_Test
View(Prostate_Low_Grade_Normal_Final@meta.data)

#Running PCA
TEST <- Prostate_Low_Grade_Normal_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 1.0)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

Prostate_Low_Grade_Normal_Clean <- TEST

DimPlot(Prostate_Low_Grade_Normal_Clean, reduction = "umap", label = TRUE)
VlnPlot(Prostate_Low_Grade_Normal_Clean, features = Tcellgenes)

Idents(Prostate_Low_Grade_Normal_Final) <- 'seurat_clusters'
QuestionCluster <- subset(Prostate_Low_Grade_Normal_Final, idents = "3")

Idents(Prostate_Low_Grade_Normal_Final) <- 'seurat_clusters'
QuestionCluster <- subset(Prostate_Low_Grade_Normal_Final, idents = "3")
cluster5=FindMarkers(Prostate_Low_Grade_Normal_Final, ident.1 = 5, min.pct = 0.25)

Prostate_Low_Grade_Normal_Clean <- subset(Prostate_Low_Grade_Normal_Clean, idents = c(0,1,2,4,5,6,7,8))
(Prostate_Low_Grade_Normal_Clean@meta.data)

##Prostate Low Grade Cancer
#Check the metadata to see what is present and what you need to add
View(Prostate_High_Grade_Tumor@meta.data)

#Editing to have common column names 
colnames(Prostate_Low_Grade_Tumor@meta.data)
colnames(Prostate_Low_Grade_Tumor@meta.data)[1] <- 'PatientID'
colnames(Prostate_Low_Grade_Tumor@meta.data)[2] <- 'nCount_originalexp'
colnames(Prostate_Low_Grade_Tumor@meta.data)[3] <- 'nFeature_originalexp'
Prostate_Low_Grade_Tumor@meta.data$PatientID %>% table()

#Changing patient names 
Idents(Prostate_Low_Grade_Tumor) <- 'PatientID'
Prostate_Low_Grade_Tumor$PatientID = Idents(Prostate_Low_Grade_Tumor)
Prostate_Low_Grade_Tumor = RenameIdents(Prostate_Low_Grade_Tumor, "SCG.PCA12.T.LG" = "Prostate_Low_Grade_002")
Prostate_Low_Grade_Tumor$PatientID = Idents(Prostate_Low_Grade_Tumor)
Prostate_Low_Grade_Tumor = RenameIdents(Prostate_Low_Grade_Tumor, "SCG.PCA17.T.LG" = "Prostate_Low_Grade_003")
Prostate_Low_Grade_Tumor$PatientID = Idents(Prostate_Low_Grade_Tumor)
Prostate_Low_Grade_Tumor = RenameIdents(Prostate_Low_Grade_Tumor, "SCG.PCA18.T.LG" = "Prostate_Low_Grade_004")
Prostate_Low_Grade_Tumor$PatientID = Idents(Prostate_Low_Grade_Tumor)
Prostate_Low_Grade_Tumor = RenameIdents(Prostate_Low_Grade_Tumor, "SCG.PCA20.T.LG" = "Prostate_Low_Grade_005")
Prostate_Low_Grade_Tumor$PatientID = Idents(Prostate_Low_Grade_Tumor)
Prostate_Low_Grade_Tumor = RenameIdents(Prostate_Low_Grade_Tumor, "SCG.PCA21.T.LG" = "Prostate_Low_Grade_006")
Prostate_Low_Grade_Tumor$PatientID = Idents(Prostate_Low_Grade_Tumor)
Prostate_Low_Grade_Tumor = RenameIdents(Prostate_Low_Grade_Tumor, "SCG.PCA3.T.LG" = "Prostate_Low_Grade_007")
Prostate_Low_Grade_Tumor$PatientID = Idents(Prostate_Low_Grade_Tumor)
Prostate_Low_Grade_Tumor = RenameIdents(Prostate_Low_Grade_Tumor, "SCG.PCA5.T.LG" = "Prostate_Low_Grade_008")
Prostate_Low_Grade_Tumor$PatientID = Idents(Prostate_Low_Grade_Tumor)
Prostate_Low_Grade_Tumor = RenameIdents(Prostate_Low_Grade_Tumor, "SCG.PCA9.T.LG" = "Prostate_Low_Grade_009")
Prostate_Low_Grade_Tumor$PatientID = Idents(Prostate_Low_Grade_Tumor)
Prostate_Low_Grade_Tumor@meta.data$PatientID %>% table()

#Adding missing columns to the metadata
Prostate_Low_Grade_Tumor_Test <- addcolumns(Prostate_Low_Grade_Tumor,"Prostate","NK","Tumor","Hirz et al., 2023 (Prostate)", "PatientID")
View(Prostate_Low_Grade_Tumor_Final@meta.data)
Prostate_Low_Grade_Tumor_Final <- Prostate_Low_Grade_Tumor_Test
View(Prostate_Low_Grade_Tumor_Final@meta.data)

#Running PCA
TEST <- Prostate_Low_Grade_Complete_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

Prostate_Low_Grade_Complete_Clean <- TEST

DimPlot(Prostate_Low_Grade_Complete_Clean, reduction = "umap", label = TRUE)

#Checking if there are a large significant population of T cells and filtering the data to remove these large populations 
VlnPlot(Prostate_Low_Grade_Complete_Clean, features = Tcellgenes)

Prostate_Low_Grade_Tumor_Clean <- Prostate_Low_Grade_Tumor_Final
(Prostate_Low_Grade_Tumor_Clean@meta.data)

###Merging normal and tumor 
Prostate_Low_Grade_Complete_Clean <- merge(x=Prostate_Low_Grade_Tumor_Clean, y=Prostate_Low_Grade_Normal_Clean)
View(Prostate_Low_Grade_Complete_Clean@meta.data)

TEST <- Prostate_Low_Grade_Complete_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

Prostate_Low_Grade_Complete_Clean <- TEST

DimPlot(Prostate_Low_Grade_Complete_Clean, reduction = "umap", label = TRUE, group.by="Condition")

#Checking if there are a large significant population of T cells and filtering the data to remove these large populations 
VlnPlot(Prostate_Low_Grade_Complete_Clean, features = Tcellgenes)

#Check how many cells per patient per condition
Prostate_Low_Grade_Complete_Clean$cellnumber <- paste(Prostate_Low_Grade_Complete_Clean$PatientID, Prostate_Low_Grade_Complete_Clean$Condition, sep = "_")
Prostate_Low_Grade_Complete_Clean@meta.data$Condition %>% table()

saveRDS(Prostate_Low_Grade_Complete_PatientsRemoved, file = "Hirz et al., 2023 (Low Grade Prostate, n=6)")
#Pancreatic Cancer ---------- 

###Reading in files
Pancreatic_Cancer_Final <- readRDS("~/Documents/areeba/Codes for analysis/PDAC_lymphoid_final.rds")
Pancreatic_Normal_Final <- readRDS("~/Documents/areeba/Codes for analysis/PDAC_Normal_Lymphoid_All.rds") 
View(Pancreatic_Normal_Final@meta.data)
View(Pancreatic_Cancer_Final@meta.data)

###Changing assay to originalexp
Pancreatic_Normal_Final[['originalexp']] = Pancreatic_Normal_Final[['RNA']]
Pancreatic_Normal_Final[['RNA']] = NULL
Pancreatic_Cancer_Final[['originalexp']] = Pancreatic_Cancer_Final[['RNA']]
Pancreatic_Cancer_Final[['RNA']] = NULL

###Subsetting in NK cells
#First check what cell types are present
Pancreatic_Normal_Final@meta.data$celltype %>% table()
Pancreatic_Cancer_Final@meta.data$celltype %>% table()
#Grab anything labeled NK except for NK-T and put into "idents"
Idents(Pancreatic_Cancer_Final) <- "celltype"
Idents(Pancreatic_Normal_Final) <- "celltype"
PancreasN = subset(Pancreatic_Normal_Final, idents='NK')
PancreasC = subset(Pancreatic_Cancer_Final, idents = c('Tissue-resident NK', 'NK'))
#Check that all cells have been subsetted
PancreasN@meta.data$celltype %>% table()
PancreasC@meta.data$celltype %>% table()

###Editing metadata
#Edit metadata to have columns we want in the right order 
#(i.e. "nCount_originalexp", "nFeature_originalexp", "cancer_type", "celltype", "Patient", "Condition","Source")
##Pancreatic Cancer
#Check the metadata to see what is present and what you need to add
View(PancreasC@meta.data)

#Editing to have common column names 
colnames(PancreasC@meta.data)
colnames(PancreasC@meta.data)[4] <- 'PatientID'
colnames(PancreasC@meta.data)[2] <- 'nCount_originalexp'
colnames(PancreasC@meta.data)[3] <- 'nFeature_originalexp'
colnames(PancreasC@meta.data)

#Changing patient names 
PancreasC@meta.data$PatientID %>% table()
Idents(PancreasC) <- 'PatientID'
PancreasC = RenameIdents(PancreasC, "Tumor1" = "Pancreas_001")
PancreasC$PatientID = Idents(PancreasC)
PancreasC = RenameIdents(PancreasC, "Tumor2" = "Pancreas_002")
PancreasC$PatientID = Idents(PancreasC)
PancreasC = RenameIdents(PancreasC, "Tumor3" = "Pancreas_003")
PancreasC$PatientID = Idents(PancreasC)
PancreasC@meta.data$PatientID %>% table()

#Adding missing columns to the metadata
PancreasC_Test <- addcolumns(PancreasC,"Pancreas","NK","Tumor","Steele et al., 2020 (PDA)", "PatientID")
View(PancreasC_Test@meta.data)
PancreasC_Final <- PancreasC_Test
View(PancreasC_Final@meta.data)

#Running PCA
TEST <- PancreasC_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 1.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

PancreasC_Final <- TEST

DimPlot(PancreasC_Final, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(PancreasC_Final, features = Tcellgenes)

Idents(PancreasC_Final) <- 'seurat_clusters'
QuestionCluster <- subset(PancreasN_Final, idents = "2")

PancreasC_Clean <- PancreasC_Final

##Pancreas Normal 
#Check the metadata to see what is present and what you need to add
View(PancreasN@meta.data)

#Editing to have common column names 
colnames(PancreasN@meta.data)
colnames(PancreasN@meta.data)[7] <- 'PatientID'
colnames(PancreasN@meta.data)[2] <- 'nCount_originalexp'
colnames(PancreasN@meta.data)[3] <- 'nFeature_originalexp'
colnames(PancreasN@meta.data)

#Changing patient names 
PancreasN@meta.data$PatientID %>% table()
Idents(PancreasN) <- 'PatientID'
PancreasN = RenameIdents(PancreasN, "Normal1" = "Pancreas_001")
PancreasN$PatientID = Idents(PancreasN)
PancreasN = RenameIdents(PancreasN, "Normal2" = "Pancreas_002")
PancreasN$PatientID = Idents(PancreasN)
PancreasN = RenameIdents(PancreasN, "Normal3" = "Pancreas_003")
PancreasN$PatientID = Idents(PancreasN)
PancreasN@meta.data$PatientID %>% table()

#Adding missing columns to the metadata
PancreasN_Test <- addcolumns(PancreasN,"Pancreas","NK","Normal","Steele et al., 2020 (PDA)", "PatientID")
View(PancreasN_Test@meta.data)
PancreasN_Test <- PancreasN_Final
View(PancreasN_Final@meta.data)

#Running PCA
TEST <- PancreasN_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 1.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

PancreasN_Final <- TEST

DimPlot(PancreasN_Final, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(PancreasN_Final, features = Tcellgenes)

Idents(PancreasN_Final) <- 'seurat_clusters'
QuestionCluster <- subset(PancreasN_Final, idents = 0)
cluster0_cleaned <- subset(QuestionCluster, idents = c(0,2,3,4,5,6,7,8,9))
PancreasN_Final_Clean <- subset(PancreasN_Final, idents =c(1,2,3))
PancreasN_Clean <- merge(x = PancreasN_Final_Clean, cluster0_cleaned)


###Merging normal and tumor 
Pancreas_Complete_Clean <- merge(x = PancreasC_Clean, y = PancreasN_Clean)
View(Pancreas_Complete_Clean@meta.data)

#Running PCA
TEST <- Pancreas_Complete_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

Pancreas_Complete_Clean <- TEST 

DimPlot(Pancreas_Complete_Clean, reduction = "umap", label = TRUE, group.by="Condition")

#Checking if there are a large significant population of T cells and filtering the data to remove these large populations 
VlnPlot(Pancreas_Complete_Clean, features = Tcellgenes)

#Check how many cells per patient per condition
Pancreas_Complete_Clean$cellnumber <- paste(Pancreas_Complete_Clean$PatientID, Pancreas_Complete_Clean$Condition, sep = "_")
Pancreas_Complete_Clean@meta.data$Condition %>% table()

saveRDS(Pancreas_Complete_PatientsRemoved, file = "Steele et al., 2020 (PDA, n = 2)")
#Liver Cancer (iCCA) ---------- 

###Reading in files
iCCA_Tumor_final <- readRDS("~/Documents/areeba/Codes for analysis/iCCA Tumor final.rds")
iCCA_normal_final <- readRDS("~/Documents/areeba/Codes for analysis/iCCA normal final.rds")

###Changing assay to originalexp
iCCA_normal_final[['originalexp']] = iCCA_normal_final[['RNA']]
iCCA_normal_final[['RNA']] = NULL
iCCA_Tumor_final[['originalexp']] = iCCA_Tumor_final[['RNA']]
iCCA_Tumor_final[['RNA']] = NULL

###Subsetting NK cells 
#First check what cell types are present
iCCA_normal_final@meta.data$celltype %>% table()
iCCA_Tumor_final@meta.data$celltype %>% table()
#Grab anything labeled NK except for NK-T and put into "idents"
Idents(iCCA_Tumor_final) <- "celltype"
Idents(iCCA_normal_final) <- "celltype"
iCCA_N = subset(iCCA_normal_final, idents= c('Circulating NK', 'Tissue-resident NK'))
iCCA_C = subset(iCCA_Tumor_final, idents = c('Tissue-resident NK', 'NK'))
#Check that all cells have been subsetted
iCCA_N@meta.data$celltype %>% table()
iCCA_C@meta.data$celltype %>% table()

###Editing the metadata
#Edit metadata to have columns we want in the right order 
#(i.e. "nCount_originalexp", "nFeature_originalexp", "cancer_type", "celltype", "Patient", "Condition","Source")
##iCCA Cancer 
#Check the metadata to see what is present and what you need to add
View(iCCA_C@meta.data)

#Editing to have common column names 
colnames(iCCA_C@meta.data)
colnames(iCCA_C@meta.data)[10] <- 'PatientID'
colnames(iCCA_C@meta.data)[2] <- 'nCount_originalexp'
colnames(iCCA_C@meta.data)[3] <- 'nFeature_originalexp'
colnames(iCCA_C@meta.data)

#Changing patient names 
iCCA_C@meta.data$PatientID %>% table()
Idents(iCCA_C) <- 'PatientID'
iCCA_C = RenameIdents(iCCA_C, "Patient 1" = "iCCA_001")
iCCA_C$PatientID = Idents(iCCA_C)
iCCA_C = RenameIdents(iCCA_C, "Patient 2" = "iCCA_002")
iCCA_C$PatientID = Idents(iCCA_C)
iCCA_C = RenameIdents(iCCA_C, "Patient 3" = "iCCA_003")
iCCA_C$PatientID = Idents(iCCA_C)
iCCA_C@meta.data$PatientID %>% table()
Idents(iCCA_C) <- 'PatientID'
iCCA_C <- subset(iCCA_C, idents=c('iCCA_001', 'iCCA_003'))

#Adding missing columns to the metadata
iCCA_C_Test <- addcolumns(iCCA_C,"Liver","NK","Tumor","Ma et al., 2022 (iCCA and HCC)", "PatientID")
View(iCCA_C_Test@meta.data)
iCCA_C_Final <- iCCA_C_Test
View(iCCA_C_Final@meta.data)

#Running PCA
TEST <- iCCA_C_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 1.0)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

iCCA_C_Clean <- TEST

DimPlot(iCCA_C_Clean, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(iCCA_C_Clean, features = Tcellgenes)


Idents(iCCA_C_Final) <- 'seurat_clusters'
cluster5 <- subset(iCCA_C_Final, idents = 5)
cluster2 <- subset(iCCA_C_Final, idents = 2)
cluster6 <- subset(iCCA_C_Final, idents = 6)
iCCA_C_Final_Clean <- subset(iCCA_C_Final, idents=c(0,1,3,4))

iCCA_C_Clean <- iCCA_C_Final_Clean

##iCCA Normal 
#Check the metadata to see what is present and what you need to add
View(iCCA_N@meta.data)

#Editing to have common column names 
colnames(iCCA_N@meta.data)
colnames(iCCA_N@meta.data)[8] <- 'PatientID'
colnames(iCCA_N@meta.data)[2] <- 'nCount_originalexp'
colnames(iCCA_N@meta.data)[3] <- 'nFeature_originalexp'
colnames(iCCA_N@meta.data)

#Changing patient names 
iCCA_N@meta.data$PatientID %>% table()
Idents(iCCA_N) <- 'PatientID'
iCCA_N = RenameIdents(iCCA_N, "Patient 1" = "iCCA_001")
iCCA_N$PatientID = Idents(iCCA_N)
iCCA_N = RenameIdents(iCCA_N, "Patient 3" = "iCCA_003")
iCCA_N$PatientID = Idents(iCCA_N)
iCCA_N@meta.data$PatientID %>% table()

#Adding missing columns to the metadata
iCCA_N_Test <- addcolumns(iCCA_N,"Liver","NK","Normal","Ma et al., 2022 (iCCA and HCC)", "PatientID")
View(iCCA_N_Test@meta.data)
iCCA_N_Final <- iCCA_N_Test
View(iCCA_N_Final@meta.data)

#Running
TEST <- iCCA_N_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 1.0)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

iCCA_N_Final <- TEST

DimPlot(iCCA_N_Final, reduction = "umap", label = TRUE)

#Checking if there are a large significant population of T cells and filtering the data to remove these large populations 
VlnPlot(iCCA_N_Final, features = Tcellgenes)

cluster1 <- subset(iCCA_N_Final, idents = 1)
cluster1_clean <- subset(cluster1, idents =c(0,2))
iCCA_N_Final_Clean <- subset(iCCA_N_Final, idents = c(0,2,3))
iCCA_N_Clean <- merge(x = iCCA_N_Final_Clean, y = cluster1_clean)
cluster3 <- subset(iCCA_N_Clean, idents = 3)

###Merging normal and tumor 
iCCA_Complete_Clean <- merge(x=iCCA_N_Clean, y=iCCA_C_Clean)
View(iCCA_Complete_Clean@meta.data)

#Running PCA
TEST <- iCCA_Complete_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:10, reduction = "pca")

iCCA_Complete_Clean <- TEST

DimPlot(iCCA_Complete_Clean, reduction = "umap", label = TRUE, group.by="Condition")

#Checking if there are a large significant population of T cells and filtering the data to remove these large populations 
VlnPlot(iCCA_Complete_Clean, features = Tcellgenes)

#Check how many cells per patient per condition
iCCA_Complete_Clean$cellnumber <- paste(iCCA_Complete_Clean$PatientID, iCCA_Complete_Clean$Condition, sep = "_")
iCCA_Complete_Clean@meta.data$cellnumber %>% table()

saveRDS(iCCA_Complete_Clean, file = "Ma et al., 2022 (iCCA, n = 2)")
#Liver Cancer (HCC) ---------- 
###Reading in files 
HCC_Tumor_Final <- readRDS("~/Documents/areeba/Codes for analysis/HCC Tumor final.rds")
HCC_Normal_Final <- readRDS("~/Documents/areeba/Codes for analysis/HCC_Normal_Final.rds")

###Changing assay to originalexp
HCC_Normal_Final[['originalexp']] = HCC_Normal_Final[['RNA']]
HCC_Normal_Final[['RNA']] = NULL
HCC_Tumor_Final[['originalexp']] = HCC_Tumor_Final[['RNA']]
HCC_Tumor_Final[['RNA']] = NULL

###Subsetting NK cells
#First check what cell types are present
HCC_Normal_Final@meta.data$celltype %>% table()
HCC_Tumor_Final@meta.data$celltype %>% table()
#Grab anything labeled NK except for NK-T and put into "idents"
Idents(HCC_Normal_Final) <- "celltype"
Idents(HCC_Tumor_Final) <- "celltype"
HCC_N = subset(HCC_Normal_Final, idents =c("Circulating NK","Tissue-resident NK"))
HCC_C = subset(HCC_Tumor_Final, idents=c("NK", "Tissue-resident NK"))
#Check that all cells have been subsetted
HCC_N@meta.data$celltype %>% table()
HCC_C@meta.data$celltype %>% table()

###Editing metadata 
##HCC normal 
#Check the metadata to see what is present and what you need to add
View(HCC_N@meta.data)

#Editing to have common column names
colnames(HCC_N@meta.data)
colnames(HCC_N@meta.data)[8] <- 'PatientID'
colnames(HCC_N@meta.data)[2] <- 'nCount_originalexp'
colnames(HCC_N@meta.data)[3] <- 'nFeature_originalexp'
colnames(HCC_N@meta.data)

#Changing patient names 
HCC_N@meta.data$PatientID %>% table()
Idents(HCC_N) <- 'PatientID'
HCC_N = RenameIdents(HCC_N, "Patient 1" = "HCC_001")
HCC_N$PatientID = Idents(HCC_N)
HCC_N = RenameIdents(HCC_N, "Patient 2" = "HCC_002")
HCC_N$PatientID = Idents(HCC_N)
HCC_N = RenameIdents(HCC_N, "Patient 3" = "HCC_003")
HCC_N$PatientID = Idents(HCC_N)
HCC_N = RenameIdents(HCC_N, "Patient 4" = "HCC_004")
HCC_N$PatientID = Idents(HCC_N)
Idents(HCC_N) <- 'PatientID'
HCC_N@meta.data$PatientID %>% table()

#Adding missing columns to the metadata
HCC_N_Test <- addcolumns(HCC_N,"Liver","NK","Normal","Ma et al., 2022 (iCCA and HCC)", "PatientID")
View(HCC_N_Test@meta.data)
HCC_N_Final <- HCC_N_Test
View(HCC_N_Final@meta.data)

#Running PCA
TEST <- HCC_N_Final

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:11)
TEST <- FindClusters(TEST, resolution = 0.5)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:11, reduction = "pca")

HCC_N_Final <- TEST

DimPlot(HCC_N_Final, reduction = "umap", label = TRUE)

#Checking if there is a large, significant population of T cells and filtering the data to remove these large populations
VlnPlot(HCC_N_Final, features = Tcellgenes)

##HCC Cancer 
#Check the metadata to see what is present and what you need to add
View(HCC_C@meta.data)

#Editing to have common column names 
colnames(HCC_C@meta.data)
colnames(HCC_C@meta.data)[8] <- 'PatientID'
colnames(HCC_C@meta.data)[2] <- 'nCount_originalexp'
colnames(HCC_C@meta.data)[3] <- 'nFeature_originalexp'
colnames(HCC_C@meta.data)

#Changing patient names 
HCC_C@meta.data$PatientID %>% table()
Idents(HCC_C) <- 'PatientID'
HCC_C = RenameIdents(HCC_C, "Patient 1" = "HCC_001")
HCC_C$PatientID = Idents(HCC_C)
HCC_C = RenameIdents(HCC_C, "Patient 2" = "HCC_002")
HCC_C$PatientID = Idents(HCC_C)
HCC_C = RenameIdents(HCC_C, "Patient 3" = "HCC_003")
HCC_C$PatientID = Idents(HCC_C)
HCC_C = RenameIdents(HCC_C, "Patient 4" = "HCC_004")
HCC_C$PatientID = Idents(HCC_C)
Idents(HCC_C) <- 'PatientID'
HCC_C@meta.data$PatientID %>% table()

#Adding missing columns to the metadata
HCC_C_Test <- addcolumns(HCC_C,"Liver","NK","Tumor","Ma et al., 2022 (iCCA and HCC)", "PatientID")
View(HCC_C_Test@meta.data)
HCC_C_Final <- HCC_C_Test
View(HCC_C_Final@meta.data)

#Running PCA
TEST <- HCC_C_Final_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:11)
TEST <- FindClusters(TEST, resolution = 1.0)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:11, reduction = "pca")

HCC_C_Final_Clean <- TEST

DimPlot(HCC_C_Final_Clean, reduction = "umap", label = TRUE)

#Checking if there are a large significant population of T cells and filtering the data to remove these large populations 
VlnPlot(HCC_C_Final_Clean, features = Tcellgenes)

Idents(HCC_C_Final) <- 'seurat_clusters'
cluster0 = subset(HCC_C_Final, idents = 0)
cluster2 = subset(HCC_C_Final, idents = 2)
cluster2_clean = subset(cluster2, idents =c(1,2))
cluster3 = subset(HCC_C_Final, idents = 3)
cluster3_clean = subset(cluster3, idents =c(1,2))
HCC_C_Final_Clean = subset(HCC_C_Final, idents=c(1,4,5,6,7,8))
HCC_C_Final_Clean <- merge(x = HCC_C_Final_Clean, y = cluster3_clean)
HCC_C_Final_Clean = subset(HCC_C_Final_Clean, idents = c(0,1,2,3,4,6))
cluster5 = subset(HCC_C_Final_Clean, idents = 5)
cluster7 = subset(HCC_C_Final_Clean, idents = 7)

###Merging normal and tumor 
HCC_Complete_Clean <- merge(x=HCC_C_Final_Clean, y=HCC_N_Final)
View(HCC_Complete_Clean@meta.data)

#Running PCA
TEST <- HCC_Complete_Clean

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:20)
TEST <- FindClusters(TEST, resolution = 1.2)
TEST <- RunUMAP(TEST, label = TRUE, dims = 1:20, reduction = "pca")

HCC_Complete_Clean <- TEST

DimPlot(HCC_Complete_PatientsRemoved, reduction = "umap", label = TRUE, group.by="Condition")

#Checking if there are a large significant population of T cells and filtering the data to remove these large populations 
VlnPlot(HCC_Complete_PatientsRemoved, features = Tcellgenes)

cluster9 = subset(HCC_Complete_Clean, idents = 9)
HCC_Complete_Clean = subset(HCC_Complete_Clean, idents =c(0,1,2,3,4,5,6,7,8))

#Check how many cells per patient per condition
HCC_Complete_Clean$cellnumber <- paste(HCC_Complete_Clean$PatientID, HCC_Complete_Clean$Condition, sep = "_")
HCC_Complete_Clean@meta.data$cellnumber %>% table()

Idents(HCC_Complete_Clean) <- "PatientID"
HCC_Complete_PatientsRemoved = subset(HCC_Complete_Clean, idents = c("HCC_001", "HCC_002", "HCC_003"))

saveRDS(HCC_Complete_Clean, file = "Ma et al., 2022 (HCC, n = 2)")
#UMAP synthesis ---------- 

###Before merging join Layers 
Lung_Complete_Clean = JoinLayers(Lung_Complete_Clean)
Kidney_Complete_Clean = JoinLayers(Kidney_Complete_Clean)
Brain_Complete_Clean = JoinLayers(Brain_Complete_Clean)
Colon_Complete_Clean = JoinLayers(Colon_Complete_Clean)
Prostate_High_Grade_Complete_Clean = JoinLayers(Prostate_High_Grade_Complete_Clean)
Prostate_Low_Grade_Complete_Clean = JoinLayers(Prostate_Low_Grade_Complete_Clean)
Pancreas_Complete_Clean = JoinLayers(Pancreas_Complete_Clean)
iCCA_Complete_Clean = JoinLayers(iCCA_Complete_Clean)
HCC_Complete_Clean = JoinLayers(HCC_Complete_Clean)
Ovarian_Complete_Clean = JoinLayers(Ovarian_Complete_Clean)
  
###Create object with just NK cells 
NK_cells = merge(x=Lung_Complete_Clean, y=c(Kidney_Complete_Clean, Brain_Complete_Clean, Colon_Complete_Clean, Prostate_High_Grade_Complete_Clean, Prostate_Low_Grade_Complete_Clean, Pancreas_Complete_Clean, iCCA_Complete_Clean, HCC_Complete_Clean, Ovarian_Complete_Clean))

#Save object 
saveRDS(NK_cells, file = "NK_cells.rds")

###Running harmony 
TEST <- NK_cells

TEST <- NormalizeData(TEST,normalization.method = "LogNormalize", scale.factor = 10000)
TEST <- FindVariableFeatures(TEST, selection.method = "vst", nfeatures = 2000)
TEST <- ScaleData(TEST)
TEST <- RunPCA(TEST,features = VariableFeatures({TEST}))
ElbowPlot(TEST)
TEST <- FindNeighbors(TEST, dims = 1:10)
TEST <- FindClusters(TEST, resolution = 1.5)
TEST <- RunHarmony(TEST, 'cancer_type')
TEST <- RunUMAP(object = TEST,  reduction = "harmony", label = TRUE, dims = 1:10)

NK_cells <- TEST

DimPlot(NK_cells, reduction = "umap", group.by = "cancer_type", repel = TRUE, pt.size = 0.5, label = FALSE)

#Checking T cell genes and NK cell genes
VlnPlot(object = NK_cells, features = Tcellgenes)
VlnPlot(object = NK_cells, features = NKcellgenes)

###Creating UMAPS 
DimPlot(NK_cells, reduction = "umap", group.by = "cancer_type", repel = TRUE, pt.size = 0.5, label = FALSE)
DimPlot(NK_cells, reduction = "umap", group.by = "Condition", repel = TRUE, pt.size = 0.5, label = FALSE)
DimPlot(NK_cells, reduction = "umap", group.by = "Source", repel = TRUE, pt.size = 0.5, label = FALSE)

#Module score ---------- 

###Using modules scores, we are able to isolate NK1, NK2 and NK3 modules from the object and analyze different features associated with them 
###Genes and modules used in the Genes used section 

##Creating module scores
NK1_module <- AddModuleScore(
  object = NK_cells,
  features = NK1_module,
  ctrl = 5,
  name = 'NK1_module'
)


NK2_module <- AddModuleScore(
  object = NK_cells,
  features = NK2_module,
  ctrl = 5,
  name = 'NK2_module'
)

NK3_module <- AddModuleScore(
  object = NK_cells,
  features = NK3_module,
  ctrl = 5,
  name = 'NK3_module'
)

ERStressmarkers_module <- AddModuleScore(
  object = NK_cells,
  features = ERstressmarkers,
  ctrl = 5,
  name = 'ERstressmarkers'
)

ERsurvivalmarkers_module <- AddModuleScore(
  object = NK_cells,
  features = ERsurvivalmarkers,
  ctrl = 5,
  name = 'ERsurvivalmarkers'
)

Cytotoxicity_module <- AddModuleScore(
  object = NK_cells,
  features = Cytotoxicity,
  ctrl = 5,
  name = 'Cytotoxicity'
)

Inflammatory_module <- AddModuleScore(
  object = NK_cells,
  features = Inflammatory,
  ctrl = 5,
  name = 'Inflammatory'
)

##Adding the module scores to NK_cells metadata
# Assuming NK1_module1 is the new column name
nk1_scores <- NK1_module@meta.data$NK1_module1
# Add the NK1_module column to the object
NK_cells@meta.data$NK1_module <- nk1_scores
# View the first few rows of the meta.data to verify
head(NK_cells@meta.data)

# Assuming N2_module1 is the new column name
nk2_scores <- NK2_module@meta.data$NK2_module1
# Add the NK2_module column to the object
NK_cells@meta.data$NK2_module <- nk2_scores
# View the first few rows of the meta.data to verify
head(NK_cells@meta.data)

# Assuming NK3_module1 is the new column name
nk3_scores <- NK3_module@meta.data$NK3_module1
# Add the NK3_module column to the object
NK_cells@meta.data$NK3_module <- nk3_scores
# View the first few rows of the meta.data to verify
head(NK_cells@meta.data)

# Assuming ERsurvivalmarkers1 is the new column name
ERsurvivalmarkers_scores <- ERsurvivalmarkers_module@meta.data$ERsurvivalmarkers1
# Add the ERsurvivalmarkers column to the object
NK_cells@meta.data$ERsurvivalmarkers <- ERsurvivalmarkers_scores
# View the first few rows of the meta.data to verify
head(NK_cells@meta.data)

# Assuming ERsurvivalmarkers1 is the new column name
ERstressmakers_scores <- ERStressmarkers_module@meta.data$ERsurvivalmarkers1
View(ERStressmarkers_module)
# Add the ERsurvivalmarkers column to the object
Paired_removed@meta.data$ERstressmarkers <- ERstressmakers_scores
# View the first few rows of the meta.data to verify
head(NK_cells@meta.data)

# Assuming Cytotoxicity1 is the new column name
Cytotoxicity_scores <- Cytotoxicity_module@meta.data$Cytotoxicity1
View(Cytotoxicity_module)
# Add the Cytotoxicity column to the object
NK_cells@meta.data$Cytotoxicity <- Cytotoxicity_scores
# View the first few rows of the meta.data to verify
head(NK_cells@meta.data)

# Assuming Inflammatory1 is the new column name
Inflammatory_scores <- Inflammatory_module@meta.data$Inflammatory1
View(Inflammatory_module)
# Add the Inflammatory column to the object
NK_cells@meta.data$Inflammatory <- Inflammatory_scores
# View the first few rows of the meta.data to verify
head(NK_cells@meta.data)

#Identites ---------- 
###Feature Plots showing the module scores
#Input the module score object and the associated column
 module_feature <- FeaturePlot(
  object = NK1,
  features = "ERstressmarkers_2",
  label = TRUE,
  pt.size = 1.5,
  order = TRUE,
  blend = FALSE
)

### Adjust the color scale for each module 
#NK1 = -1 to 2 
#NK2 = -1 to 2 
#NK3 = -1 to 1.5
module_feature + scale_color_gradient2(limits = c(-1, 2), midpoint = 0, low = "blue", mid = "white", high = "red")

### We can label each cell as NK1, NK2, or NK3 based on the intersection between where the module score was the highest and the associated Seurat clusters.
NK_cells$Identity[WhichCells(NK_cells, idents = c(20, 15, 13, 2, 11, 26, 24, 3, 8, 6, 17))] <- "NK1"
NK_cells$Identity[WhichCells(NK_cells, idents = c(16, 22, 23, 25, 12, 14, 19, 7, 9, 4, 10, 5, 18))] <- "NK2"
NK_cells$Identity[WhichCells(NK_cells, idents = c(0, 21, 1))] <- "NK3"

Idents(NK_cells) <- "Identity"
NK1 = subset(NK_cells, idents = "NK1")
NK2 = subset(NK_cells, idents = "NK2")
NK3 = subset(NK_cells, idents = "NK3")

#Change column name
colnames(NK1@meta.data)
colnames(NK1@meta.data)[34] <- 'NK1_ERstressmarkers'
#Change column name
colnames(NK2@meta.data)
colnames(NK2@meta.data)[34] <- 'NK2_ERstressmarkers'

#NK2 has no Tumor_Pancreas, so we subset the Normal out as we cannot use it in our comparasion 
NK2 <- subset(NK2, subset = Condition_cancer_type != "Normal_Pancreas")
NK2@meta.data$Condition_cancer_type %>% table()

#Analysis using labelling ---------- 
###Once each cell has been assigned an identity of NK1, NK2, or NK3, we created a Dot Plot using the NK1_module, NK2_module, and NK3_module genes to verify that the cells we have labeled are correct and expressing the expected genes.

NK1_module <- c('CTSD','CCL4','ADGRG1','CD38','CD247','CHST2','CX3CR1','KLRB1','LAIR2','IGFBP7','AKR1C3','FGFBP2','MYOM2','CLIC3','GZMB','PRF1','NKG7','SPON2')

NK2_module <- c('LTB','FOS','IL2RB','IL7R','PIK3R1','AREG','ZFP36L2','DUSP3','CD44','GPR183','CMC1','KLRC1','TCF7','TPT1','XCL2','XCL1','GZMK')

NK3_module <- c('CD3D','DSTN','ZBTB38','PPDPF','LGALS1','S100A4','IL32','VIM','PTMS','GZMH','CD3E','CCL5','KLRC2')

DotPlot(
  NK_cells,
  features = c(NK1_module, NK2_module, NK3_module),
  group.by = "label", 
  cols = "RdBu"
) + coord_flip() + RotatedAxis()

#UMAP Illustrating Distribution of NK1, NK2, and NK3
DimPlot(NK_cells, reduction = "umap", group.by = "Identity", repel = TRUE, pt.size = 0.5, label = TRUE)

#Analysis between condition and identity ---------- 
###After giving each cell an Identity, we then added more columns to the metadata to see composition of NK cells in terms of NK1, NK2 and NK3 in each tumor and make comparasions through different plots 
NK_cells$Condition_label <- paste(NK_cells$Condition, NK_cells$Identity, sep = "_")
NK_cells@meta.data$Condition_label %>% table()
NK_cells$Condition_cancer_type <- paste(NK_cells$Condition, NK_cells$cancer_type, sep = "_")
NK_cells@meta.data$Condition_cancer_type %>% table()
NK_cells$Condition_cancer_type_label <- paste(NK_cells$Condition, NK_cells$cancer_type, NK_cells$Identity, sep = "_")

###Using these new columns we created Violin Plots to see if there is a difference between Normal and Tumor of each cancer using the ER Stress genes
#Reorder the levels of each factor 
NK2$Condition_cancer_type <- factor(NK2$Condition_cancer_type, 
                                    levels = c("Normal_Brain", "Tumor_Brain","Normal_Colon", "Tumor_Colon", "Normal_Kidney", "Tumor_Kidney","Normal_Liver", "Tumor_Liver","Normal_Lung", "Tumor_Lung", "Normal_Ovarian", "Tumor_Ovarian", "Normal_Pancreas", "Tumor_Pancreas",  "Normal_Prostate", "Tumor_Prostate"))

NK1$Condition_cancer_type <- factor(NK1$Condition_cancer_type, 
                                    levels = c("Normal_Brain", "Tumor_Brain","Normal_Colon", "Tumor_Colon", "Normal_Kidney", "Tumor_Kidney","Normal_Liver", "Tumor_Liver","Normal_Lung", "Tumor_Lung", "Normal_Ovarian", "Tumor_Ovarian", "Normal_Pancreas", "Tumor_Pancreas",  "Normal_Prostate", "Tumor_Prostate"))

#List of the factors we want to compare 
my_comparisons_NK2 <- list(
  c("Normal_Brain", "Tumor_Brain"),
  c("Normal_Colon", "Tumor_Colon"),
  c("Normal_Kidney", "Tumor_Kidney"),
  c("Normal_Liver", "Tumor_Liver"),
  c("Normal_Lung", "Tumor_Lung"),
  c("Normal_Prostate", "Tumor_Prostate")
)

my_comparisons_NK1 <- list(
  c("Normal_Brain", "Tumor_Brain"),
  c("Normal_Colon", "Tumor_Colon"),
  c("Normal_Kidney", "Tumor_Kidney"),
  c("Normal_Liver", "Tumor_Liver"),
  c("Normal_Lung", "Tumor_Lung"),
  c("Normal_Ovarian", "Tumor_Ovarian"),
  c("Normal_Pancreas", "Tumor_Pancreas"),
  c("Normal_Prostate", "Tumor_Prostate")
)

#Violin Plot 
VlnPlot_scCustom(seurat_object = NK1, features = "NK1_ERstressmarkers", group.by = "Condition_cancer_type", pt.size = 0, colors_use = NULL, ggplot_default_colors = TRUE, plot_median = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means(comparisons = my_comparisons_NK1) +
  ylim(-1, 5)
#DESeq ---------- 
###pseudo-bulk workflow 
#Aggregate expression 
Aggregated_NK_cells <- AggregateExpression(NK_cells, 
                                         group.by = c("Condition", "PatientID", "cancer_type"),
                                         slot = "counts",
                                         return.seurat = TRUE)

#Extract counts matrix 
counts_matrix <- Aggregated_NK_cells[["originalexp"]]$counts
#Check to make sure that all the desired columns are present 
head(Aggregated_NK_cells)

###Making DESeq object 
#Turn counts_matrix into data frame 
DF_matrix <- as.data.frame(counts_matrix)
#Check to make sure a data frame has been produced
View(DF_matrix)

#Save the counts matrix and the metadata
write.csv (DF_matrix, file = "counts_matrix.csv", row.names =T)

write.csv (Aggregated_NK_cells@meta.data, file = "metadata_.csv")

#Extract PatientID from Aggregated object 
patients <- as.list(Aggregated_NK_cells@meta.data[,3])
View(patients)

#Set conditions for DESeq
condition <- c(rep("Normal", 64), rep("Tumor", 67))
patientID <- unlist(patients)
my_conData <- as.data.frame(condition)
my_patData <- as.data.frame(patientID)
rownames(my_conData) <- colnames(data)
my_colData <- cbind(my_conData,my_patData) 
View(my_colData)

#Check that the column name line up 
all(colnames(DF_matrix_eight_cancers_no_filter) == my_colData$PatientID)

#IF THEY DO NOT LINE UP USE THIS CODE: 
counts <- read.csv("~/Documents/areeba/Codes for analysis/counts_matrix.csv")
View(counts)

head(counts)
rownames(counts) <- counts[,1]
counts <- counts[,-1]
head(counts)

all(colnames(DF_matrix_eight_cancers_no_filter) == my_colData$PatientID)

#Make your DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(DF_matrix_eight_cancers_no_filter),
                              colData = my_colData,
                              design = ~condition+patientID)
dds <- DESeq(dds)
View(dds)

###Results and Interpreting the results
#Create results object 
res <- results(dds, contrast = c("condition", "Normal", "Tumor"))[,-c(3,4)]

#Combine normalized counts with entire DE list
normalized_counts <- round(counts(dds, normalized = TRUE),3)
gene_names <- as.data.frame(res@rownames)
pattern <- str_c("Tumor", "|", "Normal")
combined_data <- as_tibble(cbind(gene_names,res, normalized_counts))

#Generate sorted lists with the indicated cutoff values
res <- res[order(res$log2FoldChange, decreasing=TRUE),]
de_genes_pvalue <- combined_data[which(combined_data$pvalue < 0.05),]
de_genes_log2f <- combined_data[which(abs(combined_data$log2FoldChange) > 0.5 & combined_data$pvalue < 0.05),]
de_genes_cpm <- combined_data[which(combined_data$avg_cpm > 2 & combined_data$padj < 0.05),]
de_genes_padj <- combined_data[which(combined_data$padj < 0.05),]

comparisons <- c("Normal", "Tumor")

#Write output to files
write.csv (de_genes_pvalue, file = paste0(comparisons[1], "_vs_", comparisons[2], "_pvalue_cutoff2_paired.csv"), row.names =F)
write.csv (de_genes_padj, file = paste0(comparisons[1], "_vs_", comparisons[2], "_padj_cutoff1_paired.csv"), row.names =F)
write.csv (de_genes_log2f, file = paste0(comparisons[1], "_vs_", comparisons[2], "_log2f_cutoff1_paired.csv"), row.names =F)
write.csv (de_genes_cpm, file = paste0(comparisons[1], "_vs_", comparisons[2], "_cpm_cutoff1_paired.csv"), row.names =F)
write.csv (combined_data, file = paste0(comparisons[1], "_vs_", comparisons[2], "_allgenes1_paired.csv"), row.names =F)
write.table (combined_data, file = paste0(comparisons[1], "_vs_", comparisons[2], "_rank1_paired.rnk"), sep = "\t", row.names = F, quote = F)


saveRDS(Paired_removed_Kidney, file = "Paired_removed_Kidney_Liver.rds")

Paired_removed <- merge(x = Paired_removed, y = c(iCCA_Complete_Clean, HCC_Complete_PatientsRemoved))
Paired_removed <- merge(x = Paired_removed, y = list(iCCA_Complete_Clean, HCC_Complete_PatientsRemoved))
print(Paired_removed)
Idents(Paired_removed) <- "cancer_type"
Colon_Paired_removed = subset(Paired_removed, idents = "Colon")

Paired_removed_Kidney@meta.data$cancer_type %>% table()
Paired_removed_Kidney <- merge(x = Paired_removed_Kidney, y = c(iCCA_Complete_Clean, HCC_Complete_PatientsRemoved))

Idents(NK_cells) <- "cancer_type"
Ovarian_Pancreas = subset(NK_cells, idents = c("Ovarian", "Pancreas"))
saveRDS(Ovarian_Pancreas, file = "Ovarian_Pancreas.rds")
