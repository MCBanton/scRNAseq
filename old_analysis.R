##################-
### Create scRNA-seq expression matrix #### 
##################-

https://satijalab.org/seurat/pbmc3k_tutorial.html
https://satijalab.org/seurat/de_vignette.html
https://hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html#expression-qc-reads
https://davetang.org/muse/2017/08/01/getting-started-seurat/
https://bioconductor.org/packages/devel/bioc/vignettes/scater/inst/doc/vignette-dataviz.html
https://www.biostars.org/p/314881/
https://bioinformatics.stackexchange.com/questions/4470/mapping-a-list-of-cells-in-seurat-featureplot
https://www.biostars.org/p/311950/
  
### Load packages
# install.packages("tidyverse", dependencies=TRUE)
# source("https://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")
# install.packages("stringi", dependencies=TRUE)
# install.packages("dplyr", dependencies=TRUE)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("scater", version = "3.8")
# install.packages("limma")
# BiocManager::install("SingleCellExperiment", version = "3.8")
# # Install hdf5 tools to install serurat package:
# # sudo apt update
# # sudo apt-get install hdf5-tools
# # sudo apt-get install libhdf5-dev
# install.packages("Seurat")
# install.packages("mclust")
# install.packages("Matrix")

library(tidyverse)
library(scater)
library(dplyr)
library(limma)
library(SingleCellExperiment)
library(Seurat)
library(mclust)
library(Matrix)

### Create fake data
# df1<- data.frame(cbind(c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10", "E11"), c(0,0,0,3,0,10,0,1,2,0, 2)))
# df1
# colnames(df1)<- c("genes", "counts")
# df1$genes<- as.character(df1$genes)
# df2<- data.frame(cbind(c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10", "E12"), c(0,0,1,5,0,16,0,0,1,0, 0)))
# df1
# colnames(df2)<- c("genes", "counts")
# df2$genes<- as.character(df2$genes)

#########################-
### Lam et al data matrix ####
#########################-

### Set working directory ####
## My computer
setwd("/Users/Matt/Documents/Other projects and helping people/Brian/180209 Single_cell_website/Data/Lam_et_al/STAR_counts/")
## Server1
setwd("/media/data/mattb/projects/Brian_scRNAseq_website/data/lam/PRJNA358404/star_output/")

### Import STAR output and concatanate to create expression matrix ####

## List contents of directory
list_of_files<- dir()
list_of_files
list_of_files<- list_of_files[grep("*_ReadsPerGene.out.tab", dir())]
list_of_files

## Import first count table:
expression_matrix<- read.delim(list_of_files[1], header= FALSE, skip=4, colClasses = c("character", "integer", rep("NULL", 2)))
dim(expression_matrix)
head(expression_matrix)

## Get Cell/experiment name and name columns in table
file_name<- gsub("_ReadsPerGene.out.tab", "", list_of_files[1])
file_name
colnames(expression_matrix)<- c("gene", file_name)
head(expression_matrix)

## For loop to loop through all files in folder and join them together to make one large expression matrix table of counts (from STAR)
for (i in 2: length(list_of_files)){
  counts<- read.delim(list_of_files[i], header= FALSE, skip=4, colClasses = c("character", "integer", rep("NULL", 2)))
  
  file_name<- gsub("_ReadsPerGene.out.tab", "", list_of_files[i])
  colnames(counts)<- c("gene", file_name)
  expression_matrix<- full_join(expression_matrix,counts,by="gene")
  
  
}

## Look at table to make sure it is ok and does not have any missing values:
dim(expression_matrix)
head(expression_matrix)
expression_matrix[1:10, 1:10]
colnames(expression_matrix)
expression_matrix[!complete.cases(expression_matrix),]
# expression_matrix[1500,2]<- NA
# expression_matrix[!complete.cases(expression_matrix),]

###Save expression matrix ####

write.csv(expression_matrix, file = "Lam_expression_matrix_raw_counts.csv", quote = FALSE)



#########################-
### Romanov et al data matrix ####
#########################-

### Set working directory ####
## Tesla01
# setwd("/media/data/mattb/data/romanov/star_output_increased_stringency")
## dplyr does not work on Brian computer so transfered data to my server and perfomed analysis from there:

## Server1
setwd("/media/data/mattb/projects/Brian_scRNAseq_website/data/romanov/star_output_increased_stringency/")




### Import STAR output and concatanate to create expression matrix ####

## List contents of directory
list_of_files<- dir()
list_of_files
list_of_files<- list_of_files[grep("*ReadsPerGene.out.tab", dir())]
list_of_files
length(list_of_files)

## Import first count table:
expression_matrix<- read.delim(list_of_files[1], header= FALSE, skip=4, colClasses = c("character", "integer", rep("NULL", 2)))
dim(expression_matrix)
head(expression_matrix)

## Get Cell/experiment name and name columns in table
file_name<- gsub("ReadsPerGene.out.tab", "", list_of_files[1])
file_name
colnames(expression_matrix)<- c("gene", file_name)
head(expression_matrix)

## For loop to loop through all files in folder and join them together to make one large expression matrix table of counts (from STAR)
for (i in 2: length(list_of_files)){
  counts<- read.delim(list_of_files[i], header= FALSE, skip=4, colClasses = c("character", "integer", rep("NULL", 2)))
  
  file_name<- gsub("ReadsPerGene.out.tab", "", list_of_files[i])
  file_name<- gsub("_trimmed.fq.gz", "", file_name)
  colnames(counts)<- c("gene", file_name)
  expression_matrix<- full_join(expression_matrix,counts,by="gene")
  
  
}

## Look at table to make sure it is ok and does not have any missing values:
dim(expression_matrix)
head(expression_matrix)
expression_matrix[1:10, 1:10]
colnames(expression_matrix)
expression_matrix[!complete.cases(expression_matrix),]
# expression_matrix[1500,2]<- NA
# expression_matrix[!complete.cases(expression_matrix),]

###Save expression matrix ####

write.csv(expression_matrix, file = "Romanov_expression_matrix_raw_counts.csv", quote = FALSE)



###########################################-

###################-
#### Expression matrix Quality Control ####
###################-

#########################-
### Lam et al data matrix ####
#########################-

### Set working directory ####

### Import expression matrix

## Server1
setwd("/media/data/mattb/projects/Brian_scRNAseq_website/R_projects/mouse_hypothalamus_scRNAseq/")

lam_expression_matrix<- read.csv("./Lam_expression_matrix_raw_counts.csv", header = TRUE, stringsAsFactors = FALSE)
lam_expression_matrix[1:10, 1:10]
rownames(lam_expression_matrix)<-lam_expression_matrix[,1]


dim(lam_expression_matrix)
lam_expression_matrix<- lam_expression_matrix[,3:ncol(lam_expression_matrix)]
colnames(lam_expression_matrix)
dim(lam_expression_matrix)
head(lam_expression_matrix)
lam_expression_matrix[1:10, 1:10]

### Import annotations
lam_annotations<- read.delim("../GSE92707_series_matrix.txt", header = TRUE, stringsAsFactors = FALSE, skip = 31)
lam_annotations[,1]
# rownames(lam_annotations)<- make.names(lam_annotations[,1], unique=TRUE)
# rownames(lam_annotations)
lam_annotations<- lam_annotations[, 2:ncol(lam_annotations)]
head(lam_annotations)
dim(lam_annotations)
lam_annotations[1:42, 2]
colnames(lam_annotations)
### Create an sc-RNA-seq scater object

## Column names for counts need to be the same as the row names for annotation file!

lam_sce_object<- SingleCellExperiment(
  assays = list(counts = as.matrix(lam_expression_matrix)), 
  colData = lam_annotations
)

lam_sce_object<- SingleCellExperiment(
  assays = list(counts = as.matrix(lam_expression_matrix))
)

### Remove genes not expressed in any cell

keep_feature <- rowSums(counts(lam_sce_object) > 0) > 0
class(keep_feature)
length(keep_feature)
sum(keep_feature)
keep_feature
lam_sce_object <- lam_sce_object[keep_feature, ]
dim(lam_annotations)
str(lam_annotations)

### Look for spike ins and mitocondrial genes?
# isSpike(umi, "ERCC") <- grepl("^ERCC-", rownames(umi))
# isSpike(umi, "MT") <- rownames(umi) %in% 
#   c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
#     "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
#     "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
#     "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
#     "ENSG00000198840")



### Calculate the quality metrics:

# umi <- calculateQCMetrics(
#   umi,
#   feature_controls = list(
#     ERCC = isSpike(umi, "ERCC"), 
#     MT = isSpike(umi, "MT")
#   )
# )

lam_sce_object <- calculateQCMetrics(
  lam_sce_object
)

### Cell QC:Library size (total number of RNA molecules detected per sample)

lam_sce_object$total_counts

hist(
  lam_sce_object$total_counts,
  breaks = 100
)
abline(v = 25000, col = "red")

### Cell QC: Total number of unique genes detected in each sample.

lam_sce_object$total_features

hist(
  lam_sce_object$total_features,
  breaks = 100
)
abline(v = 7000, col = "red")
abline(v = 11000, col = "blue")
### Cell QC: ERCCs and MTs spike ins

# plotPhenoData(
#   umi,
#   aes_string(
#     x = "total_features",
#     y = "pct_counts_MT",
#     colour = "batch"
#   )
# )


### Filter cells:

### Manual cell Filtering:
lam_sce_object$use <- (
  # sufficient features (genes)
  filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # sufficient endogenous RNA
    filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    filter_by_MT
)

table(lam_sce_object$use)

### Automatic cell Filtering:

###########-

## Set up an example SingleCellExperiment
data("sc_example_counts")
data("sc_example_cell_info")
example_sce <- SingleCellExperiment(
  assays = list(counts = sc_example_counts), 
  colData = sc_example_cell_info
)
example_sce <- scater::normalize(example_sce)

## Examples plotting PC1 and PC2
plotPCA(example_sce)
example_sce <- plotPCA(example_sce, pca_data_input = "pdata", detect_outliers = TRUE, return_SCESet = TRUE)
example_sce$outlier

###########-
assay(lam_sce_object, "logcounts_raw") <- log2(counts(lam_sce_object) + 1)

lam_sce_object <- normalize(lam_sce_object)

plotPCA(lam_sce_object)

lam_sce_object <- plotPCA(
  lam_sce_object,
  exprs_values = "norm_counts",
  #size_by = "total_features", 
  #shape_by = "use",
  pca_data_input = "pdata",
  detect_outliers = TRUE,
  return_SCESet = TRUE
)

lam_sce_object$outlier
table(lam_sce_object$outlier)

### Venn Diagrams to compare filtering methods


### Exclude genes with technical artefacts 

plotQC(lam_sce_object, type = "highest-expression")

rownames(lam_expression_matrix)


#########################-
### Seurat ####
#########################-

#########################-
### Lam et al data Seurat ####
#########################-

### Set working directory
setwd("/Users/Matt/Documents/Other projects and helping people/Brian/180209 Single_cell_website/Data/Lam_et_al/STAR_counts/")
setwd("/media/data/mattb/projects/Brian_scRNAseq_website/R_projects/mouse_hypothalamus_scRNAseq/")

## Load expression matrix containing raw read counts
lam_expression_matrix<- read.csv("./Lam_expression_matrix_raw_counts.csv", header = TRUE, stringsAsFactors = FALSE)

## Inspect expression matrix 
dim(lam_expression_matrix)
#[1] 53801   165
#165 cells(actually 163) x 53801 genes (most probably have no expression)

lam_expression_matrix[1:5, 1:5]

## Make row names for matrix the gene names in column 2 and delete columns 1 and 2.
rownames(lam_expression_matrix)<- lam_expression_matrix[,2]
rownames(lam_expression_matrix)
ncol(lam_expression_matrix)
lam_expression_matrix<- lam_expression_matrix[,3:ncol(lam_expression_matrix)]
dim(lam_expression_matrix)
lam_expression_matrix[1:5, 1:5]

# summary of total expression per single cell
summary(colSums(lam_expression_matrix))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1436707  5550876  7622997  9307741 10841520 40076345
## Lots of reads per cell!!

# check how many genes have at least one transcript in each cell
at_least_one <- apply(lam_expression_matrix, 2, function(x) sum(x>0))
at_least_one
median(at_least_one)
# summary of genes per cell
summary(at_least_one)
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")
# The median number of genes detected per cell is 15584.

hist(colSums(lam_expression_matrix),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")
# The median total counts per cell is 7,622,997.

# manually check the number of genes detected in three or more cells
# a lot of genes are not detected in 3 or more cells
tmp <- apply(lam_expression_matrix, 1, function(x) sum(x>0))
table(tmp>=3)
# FALSE  TRUE 
# 16329 37472 
table(tmp>=2)
# FALSE  TRUE 
# 14118 39683

# Check how many cells have at least 200 detected genes?
keep <- tmp>=3
tmp <- lam_expression_matrix[keep,]
?apply # Margin 2 means columns

at_least_one <- apply(tmp, 2, function(x) sum(x>0)) #Sum number of rows with values greater than 0.

summary(at_least_one)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10882   14178   15546   15456   16688   20818
dim(tmp)
#[1] 37472   163

## Create a Seurat object to contain scRNA seq data:
seurat_lam <- CreateSeuratObject(raw.data = lam_expression_matrix, 
                                 min.cells = 3, 
                                 min.genes = 500, 
                                 project = "LAM")

seurat_lam
class(seurat_lam)
?slotNames
slotNames(seurat_lam)
seurat_lam@raw.data
seurat_lam@data
seurat_lam@scale.data
str(seurat_lam)

### QC and selecting cells for further analysis ####

## Get mouse MT genes IF for my analysis? Go to human ensembl entry and find ortholog? or go to ensembl and search for MT

### Create list of potential mitocondrial genes:

## On Ensembl website:
# Click mouse genome
# Search for "MT-*"
# Select genes on the left hand filter
# Select Mouse CL54BL6 filter (542 entries)
# Select table view filter
# Click on top right hand side button to download up to 1000 entries as table


## Import list of mouse mitocondrial genes found on Ensembl
mito_genes_table<- read.csv("mouse_MT_genes.csv", header = TRUE, stringsAsFactors = FALSE)

## Explore table
dim(mito_genes_table)
head(mito_genes_table)
mito_genes_table<- dplyr::arrange(mito_genes_table, description)
head(mito_genes_table)
## Identify mitochondrially encoded genes:
mito_genes_table<- dplyr::filter(mito_genes_table, grepl('mitochondrially encoded', description))

dim(mito_genes_table)
head(mito_genes_table)
write.csv(mito_genes_table, file = "mito_genes_table.csv")

## Make list of mitocondrial gene IDS:
mito_genes <- mito_genes_table[,1]
mito_genes
length(mito_genes)

## Identify which IDs are present in my dataset:

test<-seurat_lam@raw.data[mito_genes, ]
test
dim(test)
rownames(test)

!grepl("NA", rownames(test))
mito_genes_present<- mito_genes[!grepl("NA", rownames(test))]
dim(mito_genes_present)
mito_genes_present
seurat_lam@raw.data[mito_genes_present, ]
dim(seurat_lam@raw.data[mito_genes_present, ])

## Calculate the percentage of mitcondrial gene counts per cell
percent_mito <- Matrix::colSums(seurat_lam@raw.data[mito_genes_present, ])/Matrix::colSums(seurat_lam@raw.data)
percent_mito
summary(percent_mito)

# check out the meta data
head(seurat_lam@meta.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
seurat_lam <- AddMetaData(object = seurat_lam, metadata = percent_mito, col.name = "percent_mito")

head(seurat_lam@meta.data)

# plot number of genes, UMIs, and % mitochondria
VlnPlot(object = seurat_lam, features.plot = c("nGene", "nUMI", "percent_mito"), nCol = 3)

## Comments on plot:
# With Brian's data I am not sure if I should do any filtering?
# Threee cells have high mitochondrial percentage which may indicate contaminating cytoplasmic RNA.
# Maybe 5 outliers with unusually high number of reads (nUMI)(cutoff 2.5e+7?)
# Maybe no or 1 outlier for number of genes per cell (cutoff 20000)

### GenePlots to spot cells that have a high percentage of mitochondrial RNA and to plot the relationship between the number of unique molecules and the number of genes captured.

par(mfrow = c(1, 2))
GenePlot(object = seurat_lam, gene1 = "nUMI", gene2 = "percent_mito")
GenePlot(object = seurat_lam, gene1 = "nUMI", gene2 = "nGene")

# par(mfrow = c(1, 2))
# GenePlot(object = seurat_lam, gene1 = "nUMI", gene2 = "percent_mito", pch.use = '.')
# GenePlot(object = seurat_lam, gene1 = "nUMI", gene2 = "nGene", pch.use = '.')

## Numbers above the graphs represent the correlation coefficient

## Comments on graph, looking at the plots there appears to be three cells with very high percentage of mitocondrial DNA and 5 or 6 cells with an unusually large number of reads compared to the number of genes detected (maybe a cutoff at 2.5e+07). It might also be worth removing the two cells with less than 12000 molecules?

# manual check; I already know all cells have >200 genes
table(seurat_lam@meta.data$percent_mito < 0.06)
table(seurat_lam@meta.data$nGene > 12000)
table(seurat_lam@meta.data$nGene < 20000)
table(seurat_lam@meta.data$percent_mito < 0.06 & seurat_lam@meta.data$nGene > 12000 & seurat_lam@meta.data$nGene < 20000)

### Filter cells ####

## Decided not to filter cells for now!:

# seurat_lam <- FilterCells(object = seurat_lam, 
#                     subset.names = c("nGene", "percent.mito"), 
#                     low.thresholds = c(12000, -Inf),
#                     high.thresholds = c(20000, 0.06))
# 
# seurat_lam

seurat_lam <- FilterCells(object = seurat_lam, 
                          subset.names = c("nGene", "percent_mito"), 
                          low.thresholds = c(-Inf, -Inf),
                          high.thresholds = c(Inf, Inf))

seurat_lam
# 163 cells = as expected as no filering performed

## Graphs of filtered cells
par(mfrow = c(1, 2))
GenePlot(object = seurat_lam, gene1 = "nUMI", gene2 = "percent_mito")
GenePlot(object = seurat_lam, gene1 = "nUMI", gene2 = "nGene")

### Normalise data ####

?NormalizeData

## Plot non-normlaised data to get a feel for it:
hist(colSums(seurat_lam@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")

seurat_lam <- NormalizeData(object = seurat_lam,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
# currently there is only log-normalisation available and it might not be the best one to use according to Marioni paper.

## Plot normalised data to check it is truly normal:
# seurat_lam@data[1:5, 1:5]
# colSums(seurat_lam@data)
# class(seurat_lam@data)
# class(seurat_lam@raw.data)
# matrix(seurat_lam@data)
# dim(seurat_lam@data)
# matrix(seurat_lam@data, ncol = 163)
# colSums(matrix(seurat_lam@data, ncol = 163))


hist(colSums(matrix(seurat_lam@data, ncol = 163)),
     breaks = 100,
     main = "Total expression after normalisation",
     xlab = "Sum of expression")

## had to convert seurat_lam@data to matrix for plot to work

## Comments on plot:
# Before normalisation the expression plot per cell histogram is right skewed, after normalisation it is more normal but slightly left skewed, so not perfect. Maybe filtering to remove outlier cells (described above) would be beneficial for normalisation.

### Detection of variable genes across single cells ####

## Need to explore optimal parameters for my dataset!

?FindVariableGenes()

# the variable genes slot is empty before the analysis
seurat_lam@var.genes

?FindVariableGenes
seurat_lam <- FindVariableGenes(object = seurat_lam, 
                                mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                x.low.cutoff = 0.0125, 
                                x.high.cutoff = 3,
                                y.cutoff = 0.5,
                                num.bin = 20,
                                binning.method = "equal_width",
                                do.text = FALSE)
length(seurat_lam@var.genes)

seurat_lam <- FindVariableGenes(object = seurat_lam, 
                                mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                x.low.cutoff = 0.025, 
                                x.high.cutoff = 3,
                                y.cutoff = 0.5,
                                num.bin = 20,
                                binning.method = "equal_width",
                                do.text = FALSE)
length(seurat_lam@var.genes)

seurat_lam <- FindVariableGenes(object = seurat_lam, 
                                mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                x.low.cutoff = 0.025, 
                                x.high.cutoff = 3,
                                y.cutoff = 1,
                                num.bin = 20,
                                binning.method = "equal_width",
                                do.text = FALSE)
length(seurat_lam@var.genes)

seurat_lam <- FindVariableGenes(object = seurat_lam, 
                                mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                x.low.cutoff = 0.025, 
                                x.high.cutoff = 3,
                                y.cutoff = 1,
                                num.bin = 40,
                                binning.method = "equal_width",
                                do.text = FALSE)
length(seurat_lam@var.genes)

seurat_lam <- FindVariableGenes(object = seurat_lam, 
                                mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                x.low.cutoff = 0.0125, 
                                x.high.cutoff = 3,
                                y.cutoff = 0.5,
                                num.bin = 20,
                                binning.method = "equal_frequency")
length(seurat_lam@var.genes)

seurat_lam <- FindVariableGenes(object = seurat_lam, 
                                mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                x.low.cutoff = 0.0125, 
                                x.high.cutoff = 3,
                                y.cutoff = 0.5,
                                num.bin = 20,
                                binning.method = "equal_frequency",
                                do.text = FALSE
)
length(seurat_lam@var.genes)

## Create function/for loop to test different parameters

# pbmc <- FindVariableGenes(object = pbmc, 
#                           mean.function = ExpMean, 
#                           dispersion.function = LogVMR, 
#                           x.low.cutoff = 0.0125, 
#                           x.high.cutoff = 3,
#                           y.cutoff = 0.5)

# Tutorial graph:  average expression 0 to 6, dispersion -1 to 8. Lam data: average expression 0 to 4, dispersion -4.1 to 6. Lam data has less expression and less dispersion in general? Therefore, x.high.cutoff = 3 should be ok, x.low.cutoff = 0.0125 should be ok too, y.cutoff = 0.5, maybe a smaller number will increase the number of genes. 


# vector of variable genes
seurat_lam@var.genes
head(seurat_lam@var.genes)

# number of variable genes
length(seurat_lam@var.genes)

head(seurat_lam@hvg.info)


### Scaling data to remove unwanted sources of variation ####

# slot is empty before running ScaleData()
seurat_lam@scale.data

?ScaleData

# build linear model using nUMI and percent.mito
seurat_lam <- ScaleData(object = seurat_lam,
                        vars.to.regress = c("nUMI", "percent_mito"))

class(seurat_lam@scale.data)
seurat_lam@scale.data[1:6, 1:6]



### Perform linear dimensional reduction using PCA to identify sets of genes that  will be used for clustering later####
?RunPCA
RunPCA
?SetCalcParams
seurat_lam@calc.params$RunPCA
seurat_lam@dr
seurat_lam <- RunPCA(object = seurat_lam, 
                     pc.genes = seurat_lam@var.genes,
                     do.print = TRUE,
                     pcs.print = 1:5, 
                     genes.print = 5
)
seurat_lam@calc.params
seurat_lam@dr$pca
PrintPCAParams(seurat_lam)

# The PrintPCA() function outputs a set of genes that most strongly define a set of principal components.
PrintPCA(object = seurat_lam, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

# visualise top genes associated with principal components
VizPCA(object = seurat_lam, pcs.use = 1:2)

## The PCAPlot() function plots the principal components from a PCA; cells are coloured by their identity class according to seurat_lam@ident.

PCAPlot(object = seurat_lam, dim.1 = 1, dim.2 = 2)
PCAPlot(object = seurat_lam, dim.1 = 2, dim.2 = 3)
seurat_lam@ident


# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers. The PCA was only performed on the most variable genes, which is a subset of the dataset. The ProjectPCA step scores each gene in the dataset based on their correlation with the calculated components. This is useful because there may be genes that were not detected as variable genes in the variable gene selection step, which are still strongly correlated with cellular heterogeneity.
seurat_lam <- ProjectPCA(object = seurat_lam, do.print = FALSE)

## The PCHeatmap() function produces a heatmap based on the PCA; by default the function uses the first principal component and plots 30 genes across the number of cells specified in cells.use. Setting cells.use to a number plots the “extreme” cells on both ends of the spectrum, which dramatically speeds plotting for large datasets.

# Plotting 100 most extreme cells
PCHeatmap(object = seurat_lam, pc.use = 1, cells.use = 100, do.balanced = TRUE, label.columns = FALSE)

# Plotting 500 most extreme cells
PCHeatmap(object = seurat_lam, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

# Plotting all cells:
PCHeatmap(object = seurat_lam, pc.use = 1, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = seurat_lam, pc.use = 2, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = seurat_lam,
          pc.use = 1:12,
          #cells.use = 500,
          do.balanced = TRUE, 
          label.columns = FALSE,
          use.full = FALSE)

### Determine statistically significant principal components ####

?JackStraw

## Perform JackStraw analysis to emperically identify significant PCs
seurat_lam <- JackStraw(object = seurat_lam, num.replicate = 100, display.progress = TRUE)

# warnings()

## Visualise JackStraw analysis, dots should be above dashed line for PC to be significant
JackStrawPlot(object = seurat_lam, PCs = 1:13)
JackStrawPlot(object = seurat_lam, PCs = 1:18)
JackStrawPlot(object = seurat_lam, PCs = 1:20)
## It appears that principal components 1 to 10 are highly significant and 1 to 16 are fairly significant, while 1:20 are still significant.


##Another approach for deciding how many principal components to use is to examine the standard deviations of the principle components, which is performed by the PCElbowPlot() function. A cutoff can be drawn where there is a clear elbow in the graph.

?PCElbowPlot
PCElbowPlot(object = seurat_lam)
# It looks like an elbow would fall around principle component 16?




### Clustering cells ####

###resolution = 0.6 dims.use = 1:16 ###########
#########################-
?FindClusters
?BuildSNN
seurat_romanov@snn
seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:10, 
                               resolution = c(0.6),
                               print.output = 0,
                               save.SNN = TRUE,
                               force.recalc = TRUE)


PrintFindClustersParams(seurat_romanov)

seurat_romanov@dr
# seurat_romanov@kmeans
# seurat_romanov@snn
# seurat_romanov@spatial
# GetGeneLoadings(seurat_romanov)
# PrintDim(seurat_romanov)
# PrintPCA(seurat_romanov)

RidgePlot(object = seurat_romanov, features.plot = 'PC1')
# PrintSNNParams(object = test)
# ValidateSpecificClusters(test)
# test <- ValidateSpecificClusters(test, cluster1 = 0,
#                                        cluster2 = 2,  pc.use = 1:10)
# ValidateSpecificClusters
# ??RunClassifier
# test <- RunTSNE(object = test,
#                       dims.use = 1:16,
#                       do.fast = TRUE)
# 
# TSNEPlot(object = test, do.label = TRUE)
# test <- ValidateClusters(test, pc.use = 1:10)
# AverageDetectionRate(test, thresh.min = 0)



# seurat_lam@snn
# test@snn
# class(test@snn)
# test@snn[1:5, 1:5]
# dim_matrix<-as.matrix(test@snn)
# dim(dim_matrix)
# class(dim_matrix)
# dim_matrix[1:5,1:5]
# test@ident
# sil<- silhouette(as.integer(test@ident),test@snn)
# sil
# plot(sil)
# dim_matrix<- (1-dim_matrix)
# dim_matrix
# sil<- silhouette(as.integer(test@ident),dim_matrix)
# sil
# plot(sil)
# seurat_lam@data[1:5,1:5]
# expression_matrix<- seurat_lam@data[seurat_lam@var.genes,]
# dim(expression_matrix)
# dim_matrix2<-dist(t(expression_matrix))
# class(dim_matrix2)
# dim(dim_matrix2)
# dim_matrix2
# sil<- silhouette(as.integer(test@ident),dim_matrix2)
# sil
# plot(sil)


# library(cluster)
# scaledata<- t(expression_matrix)
# number_of_clusters <- rep(0, 20)
#repeat k-means for 1:20 and extract silhouette:
# for(i in 2:20){
#   k1to20 <- kmeans(scaledata, centers = i, nstart = 25, iter.max = 20)
#   ss <- silhouette(k1to20$cluster, dist(scaledata))
#   number_of_clusters[i] <- mean(ss[, 3])
# }
# number_of_clusters
# plot(ss)

# # Plot the  average silhouette width
# plot(1:20, number_of_clusters, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
# abline(v = which.max(number_of_clusters), lty = 2)


# expression_matrix
# ?scale
# scaledata<- t(expression_matrix)
# scaledata <- t(scale((expression_matrix)))
# rownames(scaledata)

###########################################-

# seurat_romanov <- FindClusters(object = seurat_romanov,
#                      reduction.type = "pca",
#                      dims.use = 1:10, 
#                      resolution = 1.5,
#                      print.output = 0,
#                      save.SNN = TRUE,
#                      force.recalc = TRUE)

PrintFindClustersParams(object = seurat_romanov)
seurat_romanov@snn
dim(seurat_romanov@snn)
seurat_romanov@ident
class(seurat_romanov@ident)
as.integer(seurat_romanov@ident)
# ?silhouette
# seurat_lam@snn
# as.matrix(seurat_lam@snn)
# dist(as.matrix(seurat_lam@snn))
# sil<- silhouette(as.integer(seurat_lam@ident),dist(seurat_lam@snn))
# sil<- silhouette(as.integer(test2@ident),test2@snn)
# sil
# plot(sil)
# 
# df<-as.matrix(test2@snn)
# dim(df)
# df<-(1-df)
# sil<- silhouette(as.integer(test2@ident),df)
# sil
# plot(sil)
# sil<- silhouette(as.integer(test2@ident),dist(df))
# sil
# plot(sil)
# 
# ?get_dist
# ?dist
# dist(seurat_romanov@)
### Run Non-linear dimensional reduction (tSNE) ####


seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:10,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov)
TSNEPlot(object = seurat_romanov, do.label = TRUE)


###resolution = 1.5 dims.use = 1:10 ###########
seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:10, 
                               resolution = 1.5,
                               print.output = 0,
                               save.SNN = TRUE,
                               force.recalc = TRUE)

PrintFindClustersParams(object = seurat_romanov)

### Run Non-linear dimensional reduction (tSNE) ####


seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:10,
                          do.fast = TRUE)

# TSNEPlot(object = seurat_romanov)
TSNEPlot(object = seurat_romanov, do.label = TRUE)



### resolution = 1.0 dims.use = 1:10 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:10,
                               resolution = 1.0,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:10,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 2.0 dims.use = 1:10 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:10,
                               resolution = 2.0,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:10,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 0.4 dims.use = 1:10 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:10,
                               resolution = 0.4,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:10,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)



### resolution = 0.6 dims.use = 1:18 #############


### resolution = 0.6 dims.use = 1:18 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:18,
                               resolution = 0.6,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:18,
                          do.fast = TRUE)

dim.embed

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 1.0 dims.use = 1:18 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:18,
                               resolution = 1.0,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:18,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 1.5 dims.use = 1:18 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:18,
                               resolution = 1.5,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:18,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 2.0 dims.use = 1:18 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:18,
                               resolution = 2.0,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:18,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 0.4 dims.use = 1:18 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:18,
                               resolution = 0.4,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:18,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 0.4 dims.use = 1:18 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:18,
                               resolution = 0.4,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:18,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 0.6 dims.use = 1:8 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:8,
                               resolution = 0.6,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:8,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)

### resolution = 0.6 dims.use = 1:12 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:12,
                               resolution = 0.6,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:12,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 0.8 dims.use = 1:12 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:12,
                               resolution = 0.8,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:12,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)

### resolution = 1.0 dims.use = 1:12 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:12,
                               resolution = 1.0,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:12,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)



### resolution = 0.8 dims.use = 1:12 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:12,
                               resolution = 0.8,
                               print.output = 0,
                               save.SNN = TRUE,
                               force.recalc = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:12,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)




######-

saveRDS(seurat_romanov, file = "./seurat_object_romanov.rds")






### Finding differentially expressed genes (cluster biomarkers) ####

head(seurat_romanov@ident)

## find all markers of cluster 1
? FindMarkers
cluster1_markers <- FindMarkers(object = seurat_romanov,
                                ident.1 = 1,
                                min.pct = 0.25
)
## I think min.pct = 0.25 means genes expressed in at least 25% of cells!
head(cluster1_markers)

## find all markers of cluster 2
cluster2_markers <- FindMarkers(object = seurat_romanov,
                                ident.1 = 2,
                                min.pct = 0.25
)
head(cluster2_markers)

## find all markers of cluster 3
cluster3_markers <- FindMarkers(object = seurat_romanov,
                                ident.1 = 3,
                                min.pct = 0.25
)
head(cluster3_markers)


## find all markers of cluster 4
cluster4_markers <- FindMarkers(object = seurat_romanov,
                                ident.1 = 4,
                                min.pct = 0.25
)
head(cluster4_markers)


## find all markers of cluster 0
cluster0_markers <- FindMarkers(object = seurat_romanov,
                                ident.1 = 0,
                                min.pct = 0.25
)
head(cluster0_markers)


# find markers for every cluster compared to all remaining cells, report
# only the positive ones
seurat_romanov <- FindAllMarkers(object = seurat_romanov,
                                 only.pos = TRUE,
                                 min.pct = 0.25, 
                                 thresh.use = 0.25)
head(seurat_romanov)
seurat_romanov %>% group_by(cluster) %>% top_n(2, avg_logFC)
seurat_romanov_markers_top5pos<- seurat_romanov %>% group_by(cluster) %>% top_n(5, avg_logFC)
write.csv(seurat_romanov_markers_top5pos, file = "seurat_romanov_markers_top5pos.csv")
dir()

## Try MAST and deseq2 methods for ID DE genes!

## Visualise top 5 positive genes in cluster 0
## The VlnPlot() and FeaturePlot() functions can be used to visualise marker expression
VlnPlot(object = seurat_romanov, features.plot = c("ENSMUSG00000005705", "ENSMUSG00000029819", "ENSMUSG00000021091", "ENSMUSG00000028004", "ENSMUSG00000026360"))

# you can plot raw UMI counts as well
VlnPlot(object = seurat_romanov, features.plot =c("ENSMUSG00000005705", "ENSMUSG00000029819", "ENSMUSG00000021091", "ENSMUSG00000028004", "ENSMUSG00000026360"), use.raw = TRUE, y.log = TRUE)

FeaturePlot(object = seurat_romanov, 
            features.plot =c("ENSMUSG00000005705", "ENSMUSG00000029819", "ENSMUSG00000021091", "ENSMUSG00000028004", "ENSMUSG00000026360"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne"
)

FeaturePlot(object = seurat_romanov, 
            features.plot =c("ENSMUSG00000020660"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne"
)

# FeaturePlot(object = seurat_romanov,
#             features.plot = head(row.names(my_tobit), 9),
#             cols.use = c("grey", "blue"))



seurat_lam_markers

## The DoHeatmap() function creates a heatmap of genes across all cells. Below, we use plot the top 10 marker genes for the eight clusters.

head(seurat_lam_markers)
seurat_lam_markers %>%
  group_by(cluster) %>%
  top_n(10, avg_logFC) -> top10

head(top10)

DoHeatmap(object = seurat_lam,
          genes.use = top10$gene,
          #order.by.ident = TRUE,
          slim.col.label = TRUE,
          remove.key = TRUE)

## DoHeatmap generates an expression heatmap for given cells and genes. In this case, we are plotting the top 10 markers (or all markers if less than 10) for each cluster.

# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = seurat_lam, genes.use = top10$gene, remove.key = TRUE)


### Assigning cell type identity to clusters ####













#####################################################-

#########################-
### Romanov et al data Seurat ####
#########################-

### Set working directory
setwd("/Users/Matt/Documents/Other projects and helping people/Brian/180209 Single_cell_website/Data/Romanov_et_al/star_counts/")

## Load expression matrix containing raw read counts
romanov_expression_matrix<- read.csv("./Romanov_expression_matrix_raw_counts.csv", header = TRUE, stringsAsFactors = FALSE)

## Inspect expression matrix 
dim(romanov_expression_matrix)
# [1] 53801   900
#900 cells(actually 898) x 53801 genes (most probably have no expression)

romanov_expression_matrix[1:5, 1:5]

## Make row names for matrix the gene names in column 2 and delete columns 1 and 2.
rownames(romanov_expression_matrix)<- romanov_expression_matrix[,2]
rownames(romanov_expression_matrix)
ncol(romanov_expression_matrix)
romanov_expression_matrix<- romanov_expression_matrix[,3:ncol(romanov_expression_matrix)]
dim(romanov_expression_matrix)
romanov_expression_matrix[1:5, 1:5]

# summary of total expression per single cell
summary(colSums(romanov_expression_matrix))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 9594  128307  210226  281723  347758 2125106
## Lots of reads per cell!!

# check how many genes have at least one transcript in each cell
at_least_one <- apply(romanov_expression_matrix, 2, function(x) sum(x>0))
at_least_one
median(at_least_one)
# summary of genes per cell
summary(at_least_one)
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")
# The median number of genes detected per cell is 6126.5

hist(colSums(romanov_expression_matrix),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")
# The median total counts per cell is 210,226.

# manually check the number of genes detected in three or more cells
# a lot of genes are not detected in 3 or more cells
tmp <- apply(romanov_expression_matrix, 1, function(x) sum(x>0))
table(tmp>=3)
# FALSE  TRUE 
# 18735 35066 
table(tmp>=2)
# FALSE  TRUE 
# 15918 37883 

# Check how many cells have at least 200 detected genes?
keep <- tmp>=3
tmp <- romanov_expression_matrix[keep,]

at_least_one <- apply(tmp, 2, function(x) sum(x>0)) #Sum number of rows with values greater than 0.

summary(at_least_one)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1748    5006    6119    6166    7310   13153
dim(tmp)
#[1] 35066   898

## Create a Seurat object to contain scRNA seq data:
seurat_romanov <- CreateSeuratObject(raw.data = romanov_expression_matrix, 
                                     min.cells = 3, 
                                     min.genes = 500, 
                                     project = "romanov")

seurat_romanov
# An object of class seurat in project romanov 
# 35066 genes across 898 samples.

class(seurat_romanov)
slotNames(seurat_romanov)
seurat_romanov@raw.data[1:5,1:5]
# seurat_romanov@data
seurat_romanov@scale.data
str(seurat_romanov)

### QC and selecting cells for further analysis ####

## Import list of mouse mitocondrial genes found on Ensembl
mito_genes_table<- read.csv("../../Lam_et_al/STAR_counts/mouse_MT_genes.csv", header = TRUE, stringsAsFactors = FALSE)

## Explore table
dim(mito_genes_table)
head(mito_genes_table)
mito_genes_table<- dplyr::arrange(mito_genes_table, description)
head(mito_genes_table, 25)
## Identify mitochondrially encoded genes:
mito_genes_table<- dplyr::filter(mito_genes_table, grepl('mitochondrially encoded', description))

dim(mito_genes_table)
head(mito_genes_table)
# write.csv(mito_genes_table, file = "mito_genes_table.csv")

## Make list of mitocondrial gene IDS:
mito_genes <- mito_genes_table[,1]
mito_genes
length(mito_genes)

## Identify which IDs are present in my dataset:

test<-seurat_romanov@raw.data[mito_genes, ]
test
dim(test)
rownames(test)

!grepl("NA", rownames(test))
mito_genes_present<- mito_genes[!grepl("NA", rownames(test))]
mito_genes_present
seurat_romanov@raw.data[mito_genes_present, ]
dim(seurat_romanov@raw.data[mito_genes_present, ])

## Calculate the percentage of mitcondrial gene counts per cell
percent_mito <- Matrix::colSums(seurat_romanov@raw.data[mito_genes_present, ])/Matrix::colSums(seurat_romanov@raw.data)
percent_mito
summary(percent_mito)
Matrix::colSums(seurat_romanov@raw.data)/Matrix::colSums(seurat_romanov@raw.data)
# 100% = 1

# check out the meta data
head(seurat_romanov@meta.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
seurat_romanov <- AddMetaData(object = seurat_romanov, metadata = percent_mito, col.name = "percent_mito")

head(seurat_romanov@meta.data)

# plot number of genes, UMIs, and % mitochondria
VlnPlot(object = seurat_romanov, features.plot = c("nGene", "nUMI", "percent_mito"), nCol = 3)

## Comments on plot:
# Apply filtering on Romanov data or is the deposited data already filtered??
# Seven cells have high mitochondrial percentage which may indicate contaminating cytoplasmic RNA.
# Good distribution for number of genes detected per cell 
# Maybe remove a few outliers with high numbers of reads per cell. What does next plot say?



### GenePlots to spot cells that have a high percentage of mitochondrial RNA and to plot the relationship between the number of unique molecules and the number of genes captured.

par(mfrow = c(1, 2))
GenePlot(object = seurat_romanov, gene1 = "nUMI", gene2 = "percent_mito")
GenePlot(object = seurat_romanov, gene1 = "nUMI", gene2 = "nGene")
par(mfrow = c(1, 1))
# par(mfrow = c(1, 2))
# GenePlot(object = seurat_romanov, gene1 = "nUMI", gene2 = "percent_mito", pch.use = '.')
# GenePlot(object = seurat_romanov, gene1 = "nUMI", gene2 = "nGene", pch.use = '.')

## Numbers above the graphs represent the correlation coefficient

## ## Comments on graph: Looking at the plots there appears to be seven cells with higher percentage of mitocondrial DNA about 13 cells with an unusually large number of reads compared to the number of genes detected (maybe a cutoff at 1,200,000). It might also be worth removing the cells with less than 2000 molecules?

## manual check; I already know all cells have >200 genes
## Cells with more than 15% percentage of reads belonging to mitocondrially encoded genes
table(seurat_romanov@meta.data$percent_mito < 0.15)
seurat_romanov@meta.data$percent_mito < 0.15
sum(seurat_romanov@meta.data$percent_mito > 0.15)
colnames(seurat_romanov@data)[seurat_romanov@meta.data$percent_mito > 0.15]

table(seurat_romanov@meta.data$nGene > 12000)
table(seurat_romanov@meta.data$nGene < 2000)
colnames(seurat_romanov@data)[seurat_romanov@meta.data$nGene < 2000]

table(seurat_romanov@meta.data$nUMI < 20000)
colnames(seurat_romanov@data)[seurat_romanov@meta.data$nUMI < 20000]

table(seurat_romanov@meta.data$nUMI > 1200000)
colnames(seurat_romanov@data)[seurat_romanov@meta.data$nUMI > 1200000]

table(seurat_romanov@meta.data$percent_mito < 0.15 & seurat_romanov@meta.data$nGene < 12000 & seurat_romanov@meta.data$nGene > 2000 & seurat_romanov@meta.data$nUMI > 20000 & seurat_romanov@meta.data$nUMI < 1200000)

not_outlier<- seurat_romanov@meta.data$nGene < 12000 & seurat_romanov@meta.data$nGene > 2000 & seurat_romanov@meta.data$nUMI > 20000 & seurat_romanov@meta.data$nUMI < 1200000

not_outlier
!not_outlier
sum(not_outlier)
sum(!not_outlier)
outlier<- as.data.frame(!not_outlier)
rownames(outlier)<- colnames(seurat_romanov@data)
colnames(outlier)<- "potential_outliers"
## Potentially remove 34 cells out of 864?

# stash QC stats
?AddMetaData
head(seurat_romanov@meta.data)
seurat_romanov <- AddMetaData(object = seurat_romanov, metadata = outlier, col.name = "potential_outliers")

head(seurat_romanov@meta.data)


### Filter cells ####

## Decided not to filter cells for now!:

seurat_romanov <- FilterCells(object = seurat_romanov, 
                              subset.names = c("nGene", "nUMI", "percent_mito"), 
                              low.thresholds = c(-Inf, -Inf, -Inf),
                              high.thresholds = c(Inf, Inf, Inf))

seurat_romanov
# 898 cells = as expected as no filering performed

## Graphs of filtered cells
par(mfrow = c(1, 2))
GenePlot(object = seurat_romanov, gene1 = "nUMI", gene2 = "percent_mito")
GenePlot(object = seurat_romanov, gene1 = "nUMI", gene2 = "nGene")
par(mfrow = c(1, 1))

### Normalise data ####

?NormalizeData

## Plot non-normlaised data to get a feel for it:
hist(colSums(seurat_romanov@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")

seurat_romanov <- NormalizeData(object = seurat_romanov,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000)
# currently there is only log-normalisation available and it might not be the best one to use according to Marioni paper.

hist(colSums(seurat_romanov@data),
     breaks = 100,
     main = "Total expression after normalisation",
     xlab = "Sum of expression")




### Detection of variable genes across single cells ####

## Need to explore optimal parameters for my dataset!

?FindVariableGenes()

# the variable genes slot is empty before the analysis
seurat_romanov@var.genes

seurat_romanov <- FindVariableGenes(object = seurat_romanov, 
                                    mean.function = ExpMean, 
                                    dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.0125, 
                                    x.high.cutoff = 3,
                                    y.cutoff = 0.5,
                                    num.bin = 20,
                                    binning.method = "equal_width",
                                    do.text = FALSE)
length(seurat_romanov@var.genes)

seurat_romanov <- FindVariableGenes(object = seurat_romanov, 
                                    mean.function = ExpMean, 
                                    dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.025, 
                                    x.high.cutoff = 8,
                                    y.cutoff = 0.5,
                                    num.bin = 20,
                                    binning.method = "equal_width",
                                    do.text = FALSE)
length(seurat_romanov@var.genes)

seurat_romanov <- FindVariableGenes(object = seurat_romanov, 
                                    mean.function = ExpMean, 
                                    dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.025, 
                                    x.high.cutoff = 8,
                                    y.cutoff = 1,
                                    num.bin = 20,
                                    binning.method = "equal_width",
                                    do.text = FALSE)
length(seurat_romanov@var.genes)

seurat_romanov <- FindVariableGenes(object = seurat_romanov, 
                                    mean.function = ExpMean, 
                                    dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.025, 
                                    x.high.cutoff = 8,
                                    y.cutoff = 1,
                                    num.bin = 40,
                                    binning.method = "equal_width",
                                    do.text = FALSE)
length(seurat_romanov@var.genes)

seurat_romanov <- FindVariableGenes(object = seurat_romanov, 
                                    mean.function = ExpMean, 
                                    dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.0125, 
                                    x.high.cutoff = 8,
                                    y.cutoff = 0.5,
                                    num.bin = 20,
                                    binning.method = "equal_frequency")
length(seurat_romanov@var.genes)


seurat_romanov <- FindVariableGenes(object = seurat_romanov, 
                                    mean.function = ExpMean, 
                                    dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.2, 
                                    x.high.cutoff = 8,
                                    y.cutoff = 1,
                                    num.bin = 20,
                                    binning.method = "equal_width",
                                    do.text = FALSE)
length(seurat_romanov@var.genes)

seurat_romanov <- FindVariableGenes(object = seurat_romanov, 
                                    mean.function = ExpMean, 
                                    dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.05, 
                                    x.high.cutoff = 8,
                                    y.cutoff = 1,
                                    num.bin = 20,
                                    binning.method = "equal_width",
                                    do.text = FALSE)
length(seurat_romanov@var.genes)

seurat_romanov <- FindVariableGenes(object = seurat_romanov, 
                                    mean.function = ExpMean, 
                                    dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.025, 
                                    x.high.cutoff = 8,
                                    y.cutoff = 1,
                                    num.bin = 20,
                                    binning.method = "equal_width",
                                    do.text = FALSE)
length(seurat_romanov@var.genes)

## Create function/for loop to test different parameters

# pbmc <- FindVariableGenes(object = pbmc, 
#                           mean.function = ExpMean, 
#                           dispersion.function = LogVMR, 
#                           x.low.cutoff = 0.0125, 
#                           x.high.cutoff = 3,
#                           y.cutoff = 0.5)

# Tutorial graph:  average expression 0 to 6, dispersion -1 to 8. Lam data: average expression 0 to 4, dispersion -4.1 to 6. Lam data has less expression and less dispersion in general? Therefore, x.high.cutoff = 3 should be ok, x.low.cutoff = 0.0125 should be ok too, y.cutoff = 0.5, maybe a smaller number will increase the number of genes. 


# vector of variable genes
seurat_romanov@var.genes
head(seurat_romanov@var.genes)

# number of variable genes
length(seurat_romanov@var.genes)

head(seurat_romanov@hvg.info)


### Scaling data to remove unwanted sources of variation ####

# slot is empty before running ScaleData()
seurat_romanov@scale.data

?ScaleData

# build linear model using nUMI and percent.mito
seurat_romanov <- ScaleData(object = seurat_romanov,
                            vars.to.regress = c("nUMI", "percent_mito"))

class(seurat_romanov@scale.data)
seurat_romanov@scale.data[1:6, 1:6]



### Perform linear dimensional reduction using PCA to identify sets of genes that  will be used for clustering later####
?RunPCA
# RunPCA
# ?SetCalcParams
seurat_romanov@calc.params$RunPCA
seurat_romanov@dr
seurat_romanov <- RunPCA(object = seurat_romanov, 
                         pc.genes = seurat_romanov@var.genes,
                         do.print = TRUE,
                         pcs.print = 1:5, 
                         genes.print = 5
)
seurat_romanov@calc.params
seurat_romanov@dr$pca
PrintPCAParams(seurat_romanov)

# The PrintPCA() function outputs a set of genes that most strongly define a set of principal components.
PrintPCA(object = seurat_romanov, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

# visualise top genes associated with principal components
VizPCA(object = seurat_romanov, pcs.use = 1:2)

## The PCAPlot() function plots the principal components from a PCA; cells are coloured by their identity class according to seurat_romanov@ident.

PCAPlot(object = seurat_romanov, dim.1 = 1, dim.2 = 2)
PCAPlot(object = seurat_romanov, dim.1 = 2, dim.2 = 3)
seurat_romanov@ident


# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers. The PCA was only performed on the most variable genes, which is a subset of the dataset. The ProjectPCA step scores each gene in the dataset based on their correlation with the calculated components. This is useful because there may be genes that were not detected as variable genes in the variable gene selection step, which are still strongly correlated with cellular heterogeneity.
seurat_romanov <- ProjectPCA(object = seurat_romanov, do.print = FALSE)

## The PCHeatmap() function produces a heatmap based on the PCA; by default the function uses the first principal component and plots 30 genes across the number of cells specified in cells.use. Setting cells.use to a number plots the “extreme” cells on both ends of the spectrum, which dramatically speeds plotting for large datasets.

# Plotting 100 most extreme cells
PCHeatmap(object = seurat_romanov, pc.use = 1, cells.use = 100, do.balanced = TRUE, label.columns = FALSE)

# Plotting 500 most extreme cells
PCHeatmap(object = seurat_romanov, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

# Plotting all cells:
PCHeatmap(object = seurat_romanov, pc.use = 1, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = seurat_romanov, pc.use = 2, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = seurat_romanov,
          pc.use = 1:12,
          #cells.use = 500,
          do.balanced = TRUE, 
          label.columns = FALSE,
          use.full = FALSE)




## Determine statistically significant principal components ####

?JackStraw

## Perform JackStraw analysis to emperically identify significant PCs
seurat_romanov <- JackStraw(object = seurat_romanov, num.replicate = 100, display.progress = TRUE)

# warnings()

## Visualise JackStraw analysis, dots should be above dashed line for PC to be significant
JackStrawPlot(object = seurat_romanov, PCs = 1:13)
JackStrawPlot(object = seurat_romanov, PCs = 1:18)
JackStrawPlot(object = seurat_romanov, PCs = 1:20)
## It appears that principal components 1 to 10 are highly significant, while 1:20 are still significant.


##Another approach for deciding how many principal components to use is to examine the standard deviations of the principle components, which is performed by the PCElbowPlot() function. A cutoff can be drawn where there is a clear elbow in the graph.

?PCElbowPlot
PCElbowPlot(object = seurat_romanov)
# It looks like an elbow would fall around principle component 8 or 10 or 13?






### Clustering cells ####

###resolution = 0.6 dims.use = 1:16 ###########
#########################-
?FindClusters
# FindClusters
# ?BuildSNN
# seurat_romanov@snn
# test3 <- FindClusters(object = seurat_romanov,
#                       reduction.type = "pca",
#                       dims.use = 1:10, 
#                       resolution = c(1.0,1.2,1.5,1.8, 2.0, 2.5, 3.0, 5.0),
#                       print.output = 0,
#                       save.SNN = TRUE,
#                       force.recalc = TRUE)
# 
# 
# test4 <- FindClusters2(object = seurat_romanov,
#                        reduction.type = "pca",
#                        dims.use = 1:10, 
#                        resolution = c(1.0,1.2,1.5,1.8, 2.0, 2.5, 3.0, 5.0),
#                        print.output = 0,
#                        save.SNN = TRUE,
#                        force.recalc = TRUE)
# 
# test@dr
# slotNames(test)
# test@kmeans
# test@snn
# test@spatial
# test
# test2<-GetGeneLoadings(test)
# PrintDim(test)
# PrintPCA(test)
# PCAEmbed(test)
# RefinedMapping(test)
# dim(test2)
# PrintCCAParams(seurat_romanov)
# PrintDMParams(test)
# PrintCalcVarExpRatioParams(test)
# SplitDotPlotGG(test)
# RidgePlot(object = test, features.plot = 'PC1')
# PrintSNNParams(object = test)
# ValidateSpecificClusters(test)
# test <- ValidateSpecificClusters(test, cluster1 = 0,
#                                  cluster2 = 2,  pc.use = 1:10)
# ValidateSpecificClusters
# ??RunClassifier
# test <- RunTSNE(object = test,
#                 dims.use = 1:16,
#                 do.fast = TRUE)
# 
# TSNEPlot(object = test, do.label = TRUE)
# test <- ValidateClusters(test, pc.use = 1:10)
# AverageDetectionRate(test, thresh.min = 0)
# 
# 
# 
# seurat_romanov@snn
# test@snn
# class(test@snn)
# test@snn[1:5, 1:5]
# dim_matrix<-as.matrix(test@snn)
# dim(dim_matrix)
# class(dim_matrix)
# dim_matrix[1:5,1:5]
# test@ident
# sil<- silhouette(as.integer(test@ident),test@snn)
# sil
# plot(sil)
# dim_matrix<- (1-dim_matrix)
# dim_matrix
# sil<- silhouette(as.integer(test@ident),dim_matrix)
# sil
# plot(sil)
# seurat_romanov@data[1:5,1:5]
# expression_matrix<- seurat_romanov@data[seurat_romanov@var.genes,]
# dim(expression_matrix)
# dim_matrix2<-dist(t(expression_matrix))
# class(dim_matrix2)
# dim(dim_matrix2)
# dim_matrix2
# sil<- silhouette(as.integer(test@ident),dim_matrix2)
# sil
# plot(sil)
# 
# 
# library(cluster)
# scaledata<- t(expression_matrix)
# number_of_clusters <- rep(0, 20)
# #repeat k-means for 1:20 and extract silhouette:
# for(i in 2:20){
#   k1to20 <- kmeans(scaledata, centers = i, nstart = 25, iter.max = 20)
#   ss <- silhouette(k1to20$cluster, dist(scaledata))
#   number_of_clusters[i] <- mean(ss[, 3])
# }
# number_of_clusters
# plot(ss)
# 
# # Plot the  average silhouette width
# plot(1:20, number_of_clusters, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
# abline(v = which.max(number_of_clusters), lty = 2)
# 
# 
# expression_matrix
# ?scale
# scaledata<- t(expression_matrix)
# scaledata <- t(scale((expression_matrix)))
# rownames(scaledata)

###########################################-

# seurat_romanov <- FindClusters(object = seurat_romanov,
#                      reduction.type = "pca",
#                      dims.use = 1:10, 
#                      resolution = 1.5,
#                      print.output = 0,
#                      save.SNN = TRUE,
#                      force.recalc = TRUE)

# PrintFindClustersParams(object = seurat_romanov)
# seurat_romanov@snn
# dim(seurat_romanov@snn)
# seurat_romanov@ident
# class(seurat_romanov@ident)
# as.integer(seurat_romanov@ident)
# ?silhouette
# seurat_romanov@snn
# as.matrix(seurat_romanov@snn)
# dist(as.matrix(seurat_romanov@snn))
# sil<- silhouette(as.integer(seurat_romanov@ident),dist(seurat_romanov@snn))
# sil<- silhouette(as.integer(test2@ident),test2@snn)
# sil
# plot(sil)
# 
# df<-as.matrix(test2@snn)
# dim(df)
# df<-(1-df)
# sil<- silhouette(as.integer(test2@ident),df)
# sil
# plot(sil)
# sil<- silhouette(as.integer(test2@ident),dist(df))
# sil
# plot(sil)
# 
# ?get_dist
# ?dist
# dist(seurat_romanov@)

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:10, 
                               resolution = c(0.6),
                               print.output = 0,
                               save.SNN = TRUE,
                               force.recalc = TRUE)


PrintFindClustersParams(seurat_romanov)

seurat_romanov@dr
# seurat_romanov@kmeans
# seurat_romanov@snn
# seurat_romanov@spatial
# GetGeneLoadings(seurat_romanov)
# PrintDim(seurat_romanov)
# PrintPCA(seurat_romanov)

RidgePlot(object = seurat_romanov, features.plot = 'PC1')
# PrintSNNParams(object = test)
# ValidateSpecificClusters(test)
# test <- ValidateSpecificClusters(test, cluster1 = 0,
#                                        cluster2 = 2,  pc.use = 1:10)
# ValidateSpecificClusters
# ??RunClassifier
# test <- RunTSNE(object = test,
#                       dims.use = 1:16,
#                       do.fast = TRUE)
# 
# TSNEPlot(object = test, do.label = TRUE)
# test <- ValidateClusters(test, pc.use = 1:10)
# AverageDetectionRate(test, thresh.min = 0)



# seurat_lam@snn
# test@snn
# class(test@snn)
# test@snn[1:5, 1:5]
# dim_matrix<-as.matrix(test@snn)
# dim(dim_matrix)
# class(dim_matrix)
# dim_matrix[1:5,1:5]
# test@ident
# sil<- silhouette(as.integer(test@ident),test@snn)
# sil
# plot(sil)
# dim_matrix<- (1-dim_matrix)
# dim_matrix
# sil<- silhouette(as.integer(test@ident),dim_matrix)
# sil
# plot(sil)
# seurat_lam@data[1:5,1:5]
# expression_matrix<- seurat_lam@data[seurat_lam@var.genes,]
# dim(expression_matrix)
# dim_matrix2<-dist(t(expression_matrix))
# class(dim_matrix2)
# dim(dim_matrix2)
# dim_matrix2
# sil<- silhouette(as.integer(test@ident),dim_matrix2)
# sil
# plot(sil)


# library(cluster)
# scaledata<- t(expression_matrix)
# number_of_clusters <- rep(0, 20)
#repeat k-means for 1:20 and extract silhouette:
# for(i in 2:20){
#   k1to20 <- kmeans(scaledata, centers = i, nstart = 25, iter.max = 20)
#   ss <- silhouette(k1to20$cluster, dist(scaledata))
#   number_of_clusters[i] <- mean(ss[, 3])
# }
# number_of_clusters
# plot(ss)

# # Plot the  average silhouette width
# plot(1:20, number_of_clusters, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
# abline(v = which.max(number_of_clusters), lty = 2)


# expression_matrix
# ?scale
# scaledata<- t(expression_matrix)
# scaledata <- t(scale((expression_matrix)))
# rownames(scaledata)

###########################################-

# seurat_romanov <- FindClusters(object = seurat_romanov,
#                      reduction.type = "pca",
#                      dims.use = 1:10, 
#                      resolution = 1.5,
#                      print.output = 0,
#                      save.SNN = TRUE,
#                      force.recalc = TRUE)

PrintFindClustersParams(object = seurat_romanov)
seurat_romanov@snn
dim(seurat_romanov@snn)
seurat_romanov@ident
class(seurat_romanov@ident)
as.integer(seurat_romanov@ident)
# ?silhouette
# seurat_lam@snn
# as.matrix(seurat_lam@snn)
# dist(as.matrix(seurat_lam@snn))
# sil<- silhouette(as.integer(seurat_lam@ident),dist(seurat_lam@snn))
# sil<- silhouette(as.integer(test2@ident),test2@snn)
# sil
# plot(sil)
# 
# df<-as.matrix(test2@snn)
# dim(df)
# df<-(1-df)
# sil<- silhouette(as.integer(test2@ident),df)
# sil
# plot(sil)
# sil<- silhouette(as.integer(test2@ident),dist(df))
# sil
# plot(sil)
# 
# ?get_dist
# ?dist
# dist(seurat_romanov@)
### Run Non-linear dimensional reduction (tSNE) ####


seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:10,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov)
TSNEPlot(object = seurat_romanov, do.label = TRUE)


###resolution = 1.5 dims.use = 1:10 ###########
seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:10, 
                               resolution = 1.5,
                               print.output = 0,
                               save.SNN = TRUE,
                               force.recalc = TRUE)

PrintFindClustersParams(object = seurat_romanov)

### Run Non-linear dimensional reduction (tSNE) ####


seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:10,
                          do.fast = TRUE)

# TSNEPlot(object = seurat_romanov)
TSNEPlot(object = seurat_romanov, do.label = TRUE)



### resolution = 1.0 dims.use = 1:10 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:10,
                               resolution = 1.0,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:10,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 2.0 dims.use = 1:10 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:10,
                               resolution = 2.0,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:10,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 0.4 dims.use = 1:10 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:10,
                               resolution = 0.4,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:10,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)



### resolution = 0.6 dims.use = 1:18 #############


### resolution = 0.6 dims.use = 1:18 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:18,
                               resolution = 0.6,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:18,
                          do.fast = TRUE)

dim.embed

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 1.0 dims.use = 1:18 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:18,
                               resolution = 1.0,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:18,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 1.5 dims.use = 1:18 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:18,
                               resolution = 1.5,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:18,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 2.0 dims.use = 1:18 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:18,
                               resolution = 2.0,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:18,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 0.4 dims.use = 1:18 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:18,
                               resolution = 0.4,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:18,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 0.4 dims.use = 1:18 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:18,
                               resolution = 0.4,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:18,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 0.6 dims.use = 1:8 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:8,
                               resolution = 0.6,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:8,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)

### resolution = 0.6 dims.use = 1:12 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:12,
                               resolution = 0.6,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:12,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)


### resolution = 0.8 dims.use = 1:12 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:12,
                               resolution = 0.8,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:12,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)

### resolution = 1.0 dims.use = 1:12 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:12,
                               resolution = 1.0,
                               print.output = 0,
                               save.SNN = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:12,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)



### resolution = 0.8 dims.use = 1:12 #############-

seurat_romanov <- FindClusters(object = seurat_romanov,
                               reduction.type = "pca",
                               dims.use = 1:12,
                               resolution = 0.8,
                               print.output = 0,
                               save.SNN = TRUE,
                               force.recalc = TRUE)

### Run Non-linear dimensional reduction (tSNE)
seurat_romanov <- RunTSNE(object = seurat_romanov,
                          dims.use = 1:12,
                          do.fast = TRUE)

TSNEPlot(object = seurat_romanov, do.label = TRUE)




######-

# saveRDS(seurat_romanov, file = "./seurat_object_romanov.rds")






### Finding differentially expressed genes (cluster biomarkers) ####

head(seurat_romanov@ident)

## find all markers of cluster 1
? FindMarkers
cluster1_markers <- FindMarkers(object = seurat_romanov,
                                ident.1 = 1,
                                min.pct = 0.25
)
## I think min.pct = 0.25 means genes expressed in at least 25% of cells!
head(cluster1_markers)

## find all markers of cluster 2
cluster2_markers <- FindMarkers(object = seurat_romanov,
                                ident.1 = 2,
                                min.pct = 0.25
)
head(cluster2_markers)

## find all markers of cluster 3
cluster3_markers <- FindMarkers(object = seurat_romanov,
                                ident.1 = 3,
                                min.pct = 0.25
)
head(cluster3_markers)


## find all markers of cluster 4
cluster4_markers <- FindMarkers(object = seurat_romanov,
                                ident.1 = 4,
                                min.pct = 0.25
)
head(cluster4_markers)


## find all markers of cluster 0
cluster0_markers <- FindMarkers(object = seurat_romanov,
                                ident.1 = 0,
                                min.pct = 0.25
)
head(cluster0_markers)


# find markers for every cluster compared to all remaining cells, report
# only the positive ones
romanov_markers <- FindAllMarkers(object = seurat_romanov,
                                  only.pos = TRUE,
                                  min.pct = 0.25, 
                                  thresh.use = 0.25)
head(romanov_markers)
romanov_markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
seurat_romanov_markers_top5pos<- romanov_markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
seurat_romanov_markers_top5pos
write.csv(seurat_romanov_markers_top5pos, file = "seurat_romanov_markers_top5pos.csv")
dir()

### Create a table of annotated genes for the top 5 positive genes per cluster

## Import Biomart ensembl file with all mouse gene annotations.

mouse_annotations<- read.csv("../../martquery_0629110806_339.txt", header = TRUE, stringsAsFactors = FALSE)

dim(mouse_annotations)
head(mouse_annotations)

mouse_annotations %>%
  dplyr::filter(str_detect("Gene.stable.ID", "ENSMUSG00000070880"))


mouse_annotations %>%
  dplyr::filter(grepl("ENSMUSG00000070880", "Gene.stable.ID"))

mouse_annotations[grep("ENSMUSG00000070880" ,mouse_annotations$Gene.stable.ID),]

unique(mouse_annotations$Gene.stable.ID)
length(unique(mouse_annotations$Gene.stable.ID))
# mouse_annotations_unique<- distinct(mouse_annotations, Gene.stable.ID)
# dim(mouse_annotations)
# dim(mouse_annotations_unique)
# length(unique(mouse_annotations$Gene.stable.ID))
# length(unique(mouse_annotations_unique$Gene.stable.ID))

duplicated(mouse_annotations$Gene.stable.ID)
mouse_annotations<-mouse_annotations[!duplicated(mouse_annotations$Gene.stable.ID),]
dim(mouse_annotations)


## Get mouse annotations for top 5 positive genes per cluster
head(seurat_romanov_markers_top5pos)
head(mouse_annotations)
colnames(mouse_annotations)[1]<- "gene"
colnames(mouse_annotations)
seurat_romanov_markers_top5pos_annotated<- mouse_annotations %>%
  dplyr::select(Gene.name, gene, Gene.description, Gene.type) %>%
  dplyr:: inner_join(seurat_romanov_markers_top5pos, .,  by = "gene", copy = FALSE)

dim(seurat_romanov_markers_top5pos_annotated)
head(seurat_romanov_markers_top5pos_annotated)

seurat_romanov_markers_top5pos_annotated$cluster
seurat_romanov_markers_top5pos_annotated[1:10,6:10]

write.csv(seurat_romanov_markers_top5pos_annotated, file = "seurat_romanov_markers_top5pos_annotated.csv")








## Try MAST and deseq2 methods for ID DE genes!

## Visualise top 5 positive genes in cluster 0
## The VlnPlot() and FeaturePlot() functions can be used to visualise marker expression
VlnPlot(object = seurat_romanov, features.plot = c("ENSMUSG00000005705", "ENSMUSG00000029819", "ENSMUSG00000021268", "ENSMUSG00000020660", "ENSMUSG00000005268", "ENSMUSG00000038418"))

# you can plot raw UMI counts as well
# VlnPlot(object = seurat_romanov, features.plot =c("ENSMUSG00000005705", "ENSMUSG00000029819", "ENSMUSG00000021091", "ENSMUSG00000028004", "ENSMUSG00000026360"), use.raw = TRUE, y.log = TRUE)

FeaturePlot(object = seurat_romanov, 
            features.plot =c("ENSMUSG00000005705", "ENSMUSG00000029819", "ENSMUSG00000021268", "ENSMUSG00000020660", "ENSMUSG00000005268", "ENSMUSG00000038418"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne"
)

# FeaturePlot(object = seurat_romanov, 
#             features.plot =c("ENSMUSG00000020660"), 
#             cols.use = c("grey", "blue"), 
#             reduction.use = "tsne"
)

# FeaturePlot(object = seurat_romanov,
#             features.plot = head(row.names(my_tobit), 9),
#             cols.use = c("grey", "blue"))



romanov_markers

## The DoHeatmap() function creates a heatmap of genes across all cells. Below, we use plot the top 10 marker genes for the eight clusters.

head(romanov_markers)
romanov_markers %>%
  group_by(cluster) %>%
  top_n(10, avg_logFC) -> top10

head(top10)

DoHeatmap(object = seurat_romanov,
          genes.use = top10$gene,
          #order.by.ident = TRUE,
          slim.col.label = TRUE,
          remove.key = TRUE)

## DoHeatmap generates an expression heatmap for given cells and genes. In this case, we are plotting the top 10 markers (or all markers if less than 10) for each cluster.

# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
# DoHeatmap(object = seurat_romanov, genes.use = top10$gene, remove.key = TRUE)


romanov_up_down_markers <- FindAllMarkers(object = seurat_romanov,
                                          #only.pos = TRUE,
                                          min.pct = 0.25, 
                                          thresh.use = 0.25)

romanov_up_down_markers %>%
  group_by(cluster) %>%
  top_n(10, avg_logFC) -> top10_up_down

top10_up_down
class(top10_up_down)

romanov_up_down_markers %>%
  group_by(cluster) %>%
  top_n(-10, avg_logFC) %>%
  rbind(.,top10_up_down)-> top10_up_down

top10_up_down[1:20,]

top10_up_down2<- dplyr::arrange(top10_up_down, cluster, avg_logFC)
top10_up_down2[11:20,]
top10_up_down2[21:30,]

top10_up_down2 %>% 
  dplyr::filter(cluster == 1)

top10_up_down2 %>% 
  dplyr::filter(cluster == 0)

DoHeatmap(object = seurat_romanov,
          genes.use = top10_up_down2$gene,
          #order.by.ident = TRUE,
          slim.col.label = TRUE,
          remove.key = TRUE)


### Assigning cell type identity to clusters ####
