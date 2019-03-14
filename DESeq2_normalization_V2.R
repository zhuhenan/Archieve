library(RColorBrewer) 
library(ggplot2) 
library(DESeq2)
library(pheatmap)
library(readr)
library(BiocParallel)
library(apeglm)
library(IHW)
library(ggpubr)
library(gridExtra)
library(vsn)
library(dplyr)
#library(HenanRtoolbox)

#
register(MulticoreParam(20))
rm(list=ls())
setwd("~/Documents/ProjectWorkspace/DEAnalysis/")

################################################################################################
# Preparing count matrices
################################################################################################

# Grep raw count file name
batchFiles <- grep('counts.txt', list.files("RawCounts"), value = TRUE) 
batchFiles <- batchFiles[order(nchar(batchFiles), batchFiles)]
# create and view an object with file names and full paths
(batchFiles <- file.path("RawCounts", batchFiles))
# Read all tables
raw_counts <- read.delim(batchFiles[1], row.names=1, quote="")
for (file in batchFiles) {
  raw_counts <- cbind(raw_counts, read.delim(file, row.names=1, quote=""))
}
# Filter known outliers
raw_counts <- raw_counts[,!colnames(raw_counts) %in% c("GBX1007.sGluc", "LB2840", "GS017")]
# Rename specific rows
colnames(raw_counts)[colnames(raw_counts) == "GS054.sGluc.control"] <- "GS054.sGluc.control3"
colnames(raw_counts)[colnames(raw_counts) == "GBX1032.sGluc"] <- "GBX1032.1.sGluc"
colnames(raw_counts)[colnames(raw_counts) == "GS054.sGluc.control.1"] <- "GS054.sGluc.control4"
colnames(raw_counts)[colnames(raw_counts) == "GBX1054.sGluc"] <- "GBX1054.sGluc.control"
colnames(raw_counts)[colnames(raw_counts) == "GBX1104.sGluc"] <- "GBX1104.sGluc.control"
colnames(raw_counts)[colnames(raw_counts) == "GBX2025.sGluc"] <- "GBX2025.sGluc.control"
colnames(raw_counts)[colnames(raw_counts) == "GS023.sGluc"] <- "GS023.sGluc.control"
colnames(raw_counts)[colnames(raw_counts) == "GS026.sGluc"] <- "GS026.sGluc.control"
colnames(raw_counts)[colnames(raw_counts) == "GS027.sGluc"] <- "GS027.sGluc.control"
colnames(raw_counts)[colnames(raw_counts) == "GS122.non"] <- "GS122.non.control"
colnames(raw_counts)[colnames(raw_counts) == "LB3410.CD45.Cells"] <- "LB3410.CD45"
colnames(raw_counts)[colnames(raw_counts) == "LB3570.Pre.CD45"] <- "LB3570.PRE"
colnames(raw_counts)[colnames(raw_counts) == "GBX055.sGluc"] <- "GBX1055.sGluc"
colnames(raw_counts)[colnames(raw_counts) == "LB3259"] <- "GS108"
colnames(raw_counts)[colnames(raw_counts) == "LB3410.CD45"] <- "LB3410"
colnames(raw_counts) <- gsub(".sGluc","", colnames(raw_counts))
colnames(raw_counts) <- gsub(".non","", colnames(raw_counts))
ind <- which(colnames(raw_counts) == "LB3521.CD45")
colnames(raw_counts)[colnames(raw_counts) == "LB3521"] <- "LB3521.CD45"
colnames(raw_counts)[ind] <- "LB3521"
# Other known outlier
raw_counts <- raw_counts[,!colnames(raw_counts) %in% c("GS104", "GS154", "GBX1123", "GBX1133", "GBX1169", "GS081", "LB3579", "LB3388")]
# Write raw count to table
# write.table(raw_counts, "Nathanson_allBranches_rsem_genes_raw_counts.csv", sep = ",", row.names = TRUE)

# Read meta data
metadata <- read_delim("Annotation/Nathanson_pair_annotation.csv", "\t", escape_double = FALSE, trim_ws = TRUE, na = "NA")
# Only take PT, GS, XG, SDX and XDS samples
# Control group and commented samples were removed
metadata <- metadata[metadata$Sample %in% colnames(raw_counts),]
metadata <- metadata[metadata$Type   %in% c("PT", "GS", "XG", "SDX", "XDS"),]
metadata <- metadata[(metadata$Control == 0 & metadata$Notes == ""),]
# Reformat
metadata <- data.frame(lapply(metadata, factor))
# Reorder
given_order = c("PT", "GS", "XG", "SDX", "XDS")
metadata <- metadata[order(match(metadata$Type, given_order)), ]
# Extract avaliable raw counts associating with meta data
raw_counts <- raw_counts[, as.character(metadata$Sample)]
raw_counts <- round(raw_counts)

################################################################################################
# The DESeqDataSet object, sample information and the design formula
################################################################################################

# construct the DESeqDataSet object from the matrix of counts and the sample information table
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ Batch + Pair + Type)

################################################################################################
# Exploratory analysis and visualization: Preparation
################################################################################################

# Pre-filtering the dataset
# remove the rows that have no or nearly no information about the amount of gene expression
dds <- dds[rowSums(counts(dds)) > 1, ]
# rows with only zeros, and additionally many rows with only a few fragments total
#no_var_genes <- apply(raw_counts, 1, var) == 0
#raw_counts <- raw_counts[!(no_var_genes), ]
#raw_counts <- raw_counts[rowMeans(raw_counts) > quantile(rowMeans(raw_counts), 0.33),]

# differences between cell lines and treatment (the variables in the design) 
# will contribute to the expected variance-mean trend of the experiment
# For a fully unsupervised transformation, one can set blind = TRUE
# the variance stabilizing transformation (VST) for negative binomial data with a dispersion-mean trend (Anders and Huber 2010)
vsd <- vst(dds,  blind = TRUE)
# the regularized-logarithm transformation or rlog (Love, Huber, and Anders 2014)
#rld <- rlog(dds, blind = TRUE)
# simply using the log2 function (after adding 1, to avoid taking the log of zero)
# first estimate size factors to account for sequencing depth
# then specify normalized=TRUE
log2count <- log2(counts(estimateSizeFactors(dds), normalized=TRUE)[, 1:2] + 1)

# perform both transformations and compare the meanSdPlot
meanSdPlot(assay(vsd), ranks = FALSE)
# meanSdPlot(assay(rld), ranks = FALSE)
# meanSdPlot(log2count,  ranks = FALSE)
# write vsd to table
write.table(assay(vsd), file = "DESeq2_normalization_vst.csv", quote = TRUE, row.names = TRUE, sep = ",")
write.table(assay(vsd), file = "DESeq2_normalization_log2.csv", quote = TRUE, row.names = TRUE, sep = ",")

################################################################################################
# Exploratory analysis and visualization: Heatmap of sample-to-sample distances
################################################################################################

# Heatmap of sample-to-sample distances using the variance stabilizing transformation (VST)
# Calculate Sample distances
sampleDists <- dist(t(assay(vsd)))
# Reformate the sample distance matrix
sampleDistMatrix <- as.matrix(sampleDists)
# Add the heatmap row name but remove the heatmap colname
rownames(sampleDistMatrix) <- paste(vsd$Pair, vsd$Type, vsd$Sample, sep = " - " )
colnames(sampleDistMatrix) <- NULL
# Set the colour palettle
colors <- colorRampPalette(c("blue", "green" ,"red"))( 255 )
# Generate the Heatmap of sample-to-sample distances
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         fontsize_row = 6)

################################################################################################
# Principal components analysis (PCA) of the whole transcriptome
################################################################################################

# Whole transcriptom
# The Principal components analysis (PCA) using the variance stabilizing transformation (VST)
# apply plotPCA function to return the data used for plotting rather than building the plot
pcaData_whole <- plotPCA(vsd, intgroup = c("Type"), returnData = TRUE)
# Transform the contribution of variaiation to percent
percentVar    <- round(100 * attr(pcaData_whole, "percentVar"))
# Generate the PCA plot and calculating the ellipses using the multivariate normal distribution
pcaData_PT    <- pcaData_whole[pcaData_whole$Type == c("PT"), ]
pcaData_GS    <- pcaData_whole[pcaData_whole$Type == c("GS"), ]
pcaData_XG    <- pcaData_whole[pcaData_whole$Type == c("XG"), ]
pcaPlot_whole <- ggplot(pcaData_whole, aes(x = PC1, y = PC2, color = pcaData_whole$Type)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  stat_ellipse(data = pcaData_PT, aes(x = PC1, y = PC2, color = pcaData_PT$Type)) +
  stat_ellipse(data = pcaData_GS, aes(x = PC1, y = PC2, color = pcaData_GS$Type)) +
  stat_ellipse(data = pcaData_XG, aes(x = PC1, y = PC2, color = pcaData_XG$Type))
#
pcaPlot_whole
#
rm(pcaData_PT, pcaData_GS, pcaData_XG, ind, percentVar, given_order, batchFiles)

################################################################################################
# Principal components analysis (PCA) of the glioma-intrinics set
################################################################################################

# read in Verhaak annotations of glioma-intrisic genes
verhaak_intrin <- read.table("Annotation/Verhaak_2017_glioma_intrinsic.txt",sep = "\t", header = T, stringsAsFactors = F)[, 1]
# read in our list of glioma-intrinsic genes
lfc_intrin <- read.table("Annotation/p_value_lfc_glioma_intrinsic_genes.txt", sep = "\t", header = F, stringsAsFactors = F)[,1]
# build the matrix with only the consensus gene sets
intrin_mat <- assay(vsd)[intersect(verhaak_intrin, lfc_intrin), ]
# apply plotPCA function to return the data used for plotting rather than building the plot
pcaData_intrinics <- prcomp(t(intrin_mat), scale = F)
eigs <- pcaData_intrinics$sdev^2
# select the subset for ellipses
select_pc <- as.vector(metadata[which(metadata$Type == c("PT", "GS", "XG")), 1])
type_info <- metadata[, c("Sample","Type")]
rownames(type_info) <- type_info[,1]
anno_pca <- cbind(pcaData_intrinics$x, type_info)
anno_pca_PT <- anno_pca[which(anno_pca$Type == c("PT")),]
anno_pca_GS <- anno_pca[which(anno_pca$Type == c("GS")),]
anno_pca_XG <- anno_pca[which(anno_pca$Type == c("XG")),]
# Generate the PCA plot and calculating the ellipses using the multivariate t distribution
pcaPlot_intrinics <- ggplot() + 
  geom_point(aes(x = -pcaData_intrinics$x[,1], 
                 y = -pcaData_intrinics$x[,2], 
                 color = metadata$Type), size = 3) +
  labs(x = paste0("PC1: ", round(eigs[1]/sum(eigs)*100, 2), "% variance"),
       y = paste0("PC2: ", round(eigs[2]/sum(eigs)*100, 2), "% variance"),
       colors = "Sample") +
  stat_ellipse(data = anno_pca_PT, aes(x = -PC1, y = -PC2, color = anno_pca_PT$Type)) + 
  stat_ellipse(data = anno_pca_GS, aes(x = -PC1, y = -PC2, color = anno_pca_GS$Type)) +
  stat_ellipse(data = anno_pca_XG, aes(x = -PC1, y = -PC2, color = anno_pca_XG$Type))
#
pcaPlot_intrinics
#
rm(type_info, anno_pca, anno_pca_PT, anno_pca_XG, anno_pca_GS)

################################################################################################
# PCA Pair Distances
################################################################################################

# varimax rotation
pca_div <- prcomp(t(intrin_mat), scale=F)
# Scores computed via rotated loadings
max_rot <- varimax(pca_div$rotation[,1:2], normalize = F)
varimax_scores <- pca_div$x[,1:2] %*% max_rot$rotmat
# Rename columns
colnames(varimax_scores) <- c("PC1", "PC2")
# Generate a new annotation list for varimax scores
metadata_anno <- metadata[, c(1, 3, 2)]
rownames(metadata_anno) <- metadata_anno$Sample
metadata_anno <- metadata_anno[, c(2,3)]
# Attach scores matrix with annotation
varimax_scores <- cbind(varimax_scores, metadata_anno)
pair_sets <- as.list(split(varimax_scores, f = varimax_scores$Pair))
# 

# PT vs. GS PCA divergence
PT_GS <- coordinate_transform("PT", "GS", pair_sets)
#
PT_GS_plot <- ggplot(PT_GS) + 
  geom_point(aes(x = V3, y = V4)) +
  geom_segment(aes(x = V1, y = V2, xend = V3, yend = V4))

# PT vs. XG PCA divergence
PT_XG <- coordinate_transform("PT", "XG", pair_sets)
#
PT_XG_plot <- ggplot(PT_XG) + 
  geom_point(aes(x = V3, y = V4)) +
  geom_segment(aes(x = V1, y = V2, xend = V3, yend = V4))
# 
grid.arrange(PT_GS_plot, PT_XG_plot)


# GS vs. XG and GS vs. SDX PCA divergence
GS_XG  <- coordinate_transform("GS", "XG", pair_sets)
GS_SDX <- coordinate_transform("GS", "SDX", pair_sets)
GS_SDX_plot <- ggplot() + 
  geom_point(data = GS_XG, aes(x = V3, y = V4, colour = "red")) +
  geom_segment(data = GS_XG, aes(x = V1, y = V2, xend = V3, yend = V4, colour = "red")) +
  geom_point(data = GS_SDX, aes(x = V3, y = V4, colour = "blue")) +
  geom_segment(data = GS_SDX, aes(x = V1, y = V2, xend = V3, yend = V4, colour = "blue")) +
  labs(x = "PC1", y = "PC2", colour = "Comparsion")
#


# XG vs. GS and XG vs. SDX PCA divergence
XG_GS  <- coordinate_transform("XG", "GS",  pair_sets)
XG_SDX <- coordinate_transform("XG", "SDX", pair_sets)
XG_SDX_plot <- ggplot() + 
  geom_point(  data = XG_GS,  aes(x = V3, y = V4, colour = "red")) +
  geom_point(  data = XG_SDX, aes(x = V3, y = V4, colour = "blue")) +
  geom_segment(data = XG_GS,  aes(x = V1, y = V2, xend = V3, yend = V4, colour = "red")) +
  geom_segment(data = XG_SDX, aes(x = V1, y = V2, xend = V3, yend = V4, colour = "blue")) +
  labs(x = "PC1", y = "PC2", colour = "Comparsion")
#
grid.arrange(GS_SDX_plot, XG_SDX_plot)

################################################################################################
# Differential expression analysis: Running the differential expression pipeline
################################################################################################

# The differential expression analysis in DESeq2 uses a generalized linear model
# The standard differential expression analysis steps are wrapped into a single function, DESeq
# Differential analysis
seq <- DESeq(dds, parallel = TRUE)
# Save DESeq object
saveRDS(seq, file = "RDS/Nathanson_DESeq2_PT_GS_XG_SDX_XDS.rds")
# The results function contains a number of arguments to customize the results table which is generated
# The results function automatically performs independent filtering based on the mean of normalized counts for each gene
# Optimizing the number of genes which will have an adjusted p value below a given FDR cutoff, alpha
# A generalization of the idea of p value filtering is to weight hypotheses to optimize power.
# The use of the method of Independent Hypothesis Weighting (IHW) for p value adjustment of DESeq2 results

################################################################################################
# Differential expression analysis: Building the results table
################################################################################################

# Extract the estimated log2 fold changes and p values for the last variable in the design formula
# Extract the results table for leaving off the contrast argument to extract the comparison of the two levels of "dex"Type"
# Only consider results with FDR < 0.05
# Specified the coefficient or contrast we want to build a results table for PT vs. REST
pt_xg_resIHW   <- results(seq, alpha = 0.05, contrast = c("Type", "PT", "XG"),   parallel = TRUE, filterFun=ihw)
pt_gs_resIHW   <- results(seq, alpha = 0.05, contrast = c("Type", "PT", "GS"),   parallel = TRUE, filterFun=ihw)
pt_sdx_resIHW  <- results(seq, alpha = 0.05, contrast = c("Type", "PT", "SDX"),  parallel = TRUE, filterFun=ihw)
pt_xds_resIHW  <- results(seq, alpha = 0.05, contrast = c("Type", "PT", "XDS"),  parallel = TRUE, filterFun=ihw)
# Specified the coefficient or contrast we want to build a results table for XG vs. REST
xg_gs_resIHW   <- results(seq, alpha = 0.05, contrast = c("Type", "XG", "GS"),   parallel = TRUE, filterFun=ihw)
xg_sdx_resIHW  <- results(seq, alpha = 0.05, contrast = c("Type", "XG", "SDX"),  parallel = TRUE, filterFun=ihw)
xg_xds_resIHW  <- results(seq, alpha = 0.05, contrast = c("Type", "XG", "XDS"),  parallel = TRUE, filterFun=ihw)
# Specified the coefficient or contrast we want to build a results table for GS vs. REST
gs_sdx_resIHW  <- results(seq, alpha = 0.05, contrast = c("Type", "GS", "SDX"),  parallel = TRUE, filterFun=ihw)
gs_xds_resIHW  <- results(seq, alpha = 0.05, contrast = c("Type", "GS", "XDS"),  parallel = TRUE, filterFun=ihw)
# Specified the coefficient or contrast we want to build a results table for SDX vs. REST
sdx_xds_resIHW <- results(seq, alpha = 0.05, contrast = c("Type", "SDX", "XDS"), parallel = TRUE, filterFun=ihw)

# Save rds
saveRDS(pt_xg_resIHW,   file = "RDS/pt_xg_resIHW.rds"  )
saveRDS(pt_gs_resIHW,   file = "RDS/pt_gs_resIHW.rds"  )
saveRDS(pt_sdx_resIHW,  file = "RDS/pt_sdx_resIHW.rds" )
saveRDS(pt_xds_resIHW,  file = "RDS/pt_xds_resIHW.rds" )
saveRDS(xg_gs_resIHW,   file = "RDS/xg_gs_resIHW.rds"  )
saveRDS(xg_sdx_resIHW,  file = "RDS/xg_sdx_resIHW.rds" )
saveRDS(xg_xds_resIHW,  file = "RDS/xg_xds_resIHW.rds" )
saveRDS(gs_sdx_resIHW,  file = "RDS/gs_sdx_resIHW.rds" )
saveRDS(gs_xds_resIHW,  file = "RDS/gs_xds_resIHW.rds" )
saveRDS(sdx_xds_resIHW, file = "RDS/sdx_xds_resIHW.rds")

################################################################################################
# Differential expression analysis: Road rds
################################################################################################
rds_files <- list.files(path = "~/Documents/ProjectWorkspace/DEAnalysis/RDS", pattern = "*.rds", full.names = TRUE)
rds_list  <- lapply(rds_files, function(x) readRDS(x))
rds_names <- lapply(rds_files, function(x) sub('\\..*$', '', basename(x)))
names(rds_list) <- rds_names

################################################################################################
# Differential expression analysis: Exploring and exporting results
################################################################################################

# MA-Plot From Means And Log Fold Changes
# Make MA-plot which is a scatter plot of log2 fold changes (on the y-axis) versus the mean expression signal (on the x-axis)
# Colour palette: up - Red; down - Blue
# FDR < 0.05 & log2 fold change > 2
pal = c("#B31B21", "#1465AC", "darkgray")
# MA plots
pt_gs_MA  <-  ggmaplot(rds_list$pt_gs_resIHW,  main = expression("Group PT" %->% "Group GS"), 
                       fdr = 0.05, fc = 1, size = 0.4, palette = pal, 
                       genenames = as.vector(rds_list$pt_gs_resIHW$name), legend = "top", top = 0, 
                       font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                       ggtheme = ggplot2::theme_minimal())
pt_xg_MA  <-  ggmaplot(rds_list$pt_xg_resIHW,  main = expression("Group PT" %->% "Group XG"),
                       fdr = 0.05, fc = 1, size = 0.4, palette = pal, 
                       genenames = as.vector(rds_list$pt_xg_resIHW$name), legend = "top", top = 0, 
                       font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                       ggtheme = ggplot2::theme_minimal())
pt_sdx_MA <-  ggmaplot(rds_list$pt_sdx_resIHW,  main = expression("Group PT" %->% "Group SDX"),
                       fdr = 0.05, fc = 1, size = 0.4, palette = pal, 
                       genenames = as.vector(rds_list$pt_sdx_resIHW$name), legend = "top", top = 0, 
                       font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                       ggtheme = ggplot2::theme_minimal())
pt_xds_MA <-  ggmaplot(rds_list$pt_xds_resIHW, main = expression("Group PT" %->% "Group XDS"),   
                       fdr = 0.05, fc = 1, size = 0.4, palette = pal,
                       genenames = as.vector(rds_list$pt_xds_resIHW$name), legend = "top", top = 0, 
                       font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                       ggtheme = ggplot2::theme_minimal())
xg_gs_MA  <-  ggmaplot(rds_list$xg_gs_resIHW,    main = expression("Group XG" %->% "Group GS"),  
                       fdr = 0.05, fc = 1, size = 0.4, palette = pal, 
                       genenames = as.vector(rds_list$xg_gs_resIHW$name), legend = "top", top = 0, 
                       font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                       ggtheme = ggplot2::theme_minimal())
xg_sdx_MA  <- ggmaplot(rds_list$xg_sdx_resIHW,  main = expression("Group XG" %->% "Group SDX"),  
                       fdr = 0.05, fc = 1, size = 0.4, palette = pal, 
                       genenames = as.vector(rds_list$xg_sdx_resIHW$name), legend = "top", top = 0, 
                       font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                       ggtheme = ggplot2::theme_minimal())
xg_xds_MA  <- ggmaplot(rds_list$xg_xds_resIHW,  main = expression("Group XG" %->% "Group XDS"),  
                       fdr = 0.05, fc = 1, size = 0.4, palette = pal, 
                       genenames = as.vector(rds_list$xg_xds_resIHW$name), legend = "top", top = 0, 
                       font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                       ggtheme = ggplot2::theme_minimal())
gs_sdx_MA  <- ggmaplot(rds_list$gs_sdx_resIHW,  main = expression("Group GS" %->% "Group SDX"),  
                       fdr = 0.05, fc = 1, size = 0.4, palette = pal, 
                       genenames = as.vector(rds_list$gs_sdx_resIHW$name), legend = "top", top = 0, 
                       font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                       ggtheme = ggplot2::theme_minimal())
gs_xds_MA  <- ggmaplot(rds_list$gs_xds_resIHW,  main = expression("Group GS" %->% "Group XDS"),  
                       fdr = 0.05, fc = 1, size = 0.4, palette = pal, 
                       genenames = as.vector(rds_list$gs_xds_resIHW$name), legend = "top", top = 0, 
                       font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                       ggtheme = ggplot2::theme_minimal())
sdx_xds_MA <- ggmaplot(rds_list$sdx_xds_resIHW, main = expression("Group SDX" %->% "Group XDS"), 
                       fdr = 0.05, fc = 1, size = 0.4, palette = pal, 
                       genenames = as.vector(rds_list$sdx_xds_resIHW$name), legend = "top", top = 0, 
                       font.label = c("bold", 11), font.legend = "bold", font.main = "bold",
                       ggtheme = ggplot2::theme_minimal())

################################################################################################
# Differential expression analysis: Exploring and exporting results
################################################################################################

# Group all MA plot data
all_MA_data <- list(pt_gs_MA$data, pt_xg_MA$data, pt_sdx_MA$data, xg_gs_MA$data, xg_sdx_MA$data, gs_sdx_MA$data)
names(all_MA_data) <- c("pt_gs", "pt_xg", "pt_sdx", "xg_gs", "xg_sdx", "gs_sdx");
# Grep all signifcant genes form the MA plots - including both up- and down- regulated genes
all_MA_up   <- lapply(all_MA_data, function(x) levels(droplevels(x[ grep("Up",   x$sig),  ]$name)))
all_MA_down <- lapply(all_MA_data, function(x) levels(droplevels(x[ grep("Down", x$sig),  ]$name)))
# Make a list with all significant gene names
all_sig_genes <- unique(unlist(all_MA_up, use.names = FALSE), unlist(all_MA_down, use.names = FALSE))
# Build a gene list to 
gene_list <- lapply(1:length(all_sig_genes), function(x) matrix(rep(0, 16), nrow=4, ncol=4, dimnames = list(c("PT", "XG", "GS", "SDX"), c("PT", "XG", "GS", "SDX"))))
names(gene_list) <- all_sig_genes

# Assign all gene with signiture matrix with four condition
for (gene in all_sig_genes) {
  # PT vs GS
  if (gene %in% all_MA_up$pt_gs) {
    gene_list[[gene]]["PT", "GS"] <- 1
  } else if (gene %in% all_MA_down$pt_gs) {
    gene_list[[gene]]["GS", "PT"] <- 1
  }
  # PT vs XG
  if (gene %in% all_MA_up$pt_xg) {
    gene_list[[gene]]["PT", "XG"] <- 1
  } else if (gene %in% all_MA_down$pt_xg) {
    gene_list[[gene]]["XG", "PT"] <- 1
  }
  # PT vs SDX
  if (gene %in% all_MA_up$pt_sdx) {
    gene_list[[gene]]["PT", "SDX"] <- 1
  } else if (gene %in% all_MA_down$pt_sdx) {
    gene_list[[gene]]["SDX", "PT"] <- 1
  }
  # XG vs GS
  if (gene %in% all_MA_up$xg_gs) {
    gene_list[[gene]]["XG", "GS"] <- 1
  } else if (gene %in% all_MA_down$pt_sdx) {
    gene_list[[gene]]["GS", "XG"] <- 1
  }
  # XG vs SDX
  if (gene %in% all_MA_up$xg_sdx) {
    gene_list[[gene]]["XG", "SDX"] <- 1
  } else if (gene %in% all_MA_down$pt_sdx) {
    gene_list[[gene]]["SDX", "XG"] <- 1
  }
  # GS vs SDX
  if (gene %in% all_MA_up$gs_sdx) {
    gene_list[[gene]]["GS", "SDX"] <- 1
  } else if (gene %in% all_MA_down$pt_sdx) {
    gene_list[[gene]]["SDX", "GS"] <- 1
  }
}

#
load(file = "Annotation/Signiture_matrix.RData")
gene_summary <- data.frame(matrix(nrow = 0, ncol = length(names(Signiture_matrix))))
names(gene_summary) <- c(names(Signiture_matrix))

#
j <- 1
for (matrix in Signiture_matrix) {
  i <- 1
  for (gene in all_sig_genes) {
    if (identical(gene_list[[gene]], matrix)) {
      gene_summary[i, j] <- gene
      i <- i + 1
    }
  }
  j <- j + 1
}

############################################################################################ 
# PT_only, PT_low, SDX_only, GS_low, GS_only, PT_XG_High, GS_SDX_High, PT_SDX_High
#
pal <- colorRampPalette(rev(brewer.pal("RdBu", n = 11)[c(1,2,3,4,6,8,9,10,11)]))(100)
#
PT_only <- gene_summary$PT_only
heatmap_mat <- assay(vsd)[PT_only,]
pheatmap(heatmap_mat, scale = "row", color = pal, gaps_col = c(42, 64, 106, 112),
         cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T)
#
XG_only <- gene_summary$XG_only
heatmap_mat <- assay(vsd)[XG_only,]
pheatmap(heatmap_mat, scale = "row", color = pal, gaps_col = c(42, 64, 106, 112),
         cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T)
#
heatmap_mat <- assay(vsd)[SDX_only,]
pheatmap(heatmap_mat, scale = "row", color = pal, gaps_col = c(42, 64, 106, 112),
         cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T, annotation_col = sample_df, 
         annotation_names_col = T, annotation_legend = T, annotation_colors = anno_cols)
#
heatmap_mat <- assay(vsd)[GS_low,]
pheatmap(heatmap_mat, scale = "row", color = pal, gaps_col = c(42, 64, 106, 112),
         cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T, annotation_col = sample_df, 
         annotation_names_col = T, annotation_legend = T, annotation_colors = anno_cols)
#
heatmap_mat <- assay(vsd)[GS_only,]
pheatmap(heatmap_mat, scale = "row", color = pal, gaps_col = c(42, 64, 106, 112),
         cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T, annotation_col = sample_df, 
         annotation_names_col = T, annotation_legend = T, annotation_colors = anno_cols)
#
heatmap_mat <- assay(vsd)[PT_XG_High,]
pheatmap(heatmap_mat, scale = "row", color = pal, gaps_col = c(42, 64, 106, 112),
         cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T, annotation_col = sample_df, 
         annotation_names_col = T, annotation_legend = T, annotation_colors = anno_cols)
#
heatmap_mat <- assay(vsd)[PT_SDX_High,]
pheatmap(heatmap_mat, scale = "row", color = pal, gaps_col = c(42, 64, 106, 112),
         cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T, annotation_col = sample_df, 
         annotation_names_col = T, annotation_legend = T, annotation_colors = anno_cols)
#
heatmap_mat <- assay(vsd)[GS_SDX_High,]
pheatmap(heatmap_mat, scale = "row", color = pal, gaps_col = c(42, 64, 106, 112),
         cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T, annotation_col = sample_df, 
         annotation_names_col = T, annotation_legend = T, annotation_colors = anno_cols)
#
full_list <- c(PT_only, PT_SDX_High, PT_XG_High, GS_low, GS_only, PT_low, GS_SDX_High, SDX_only)
heatmap_mat <- assay(vsd)[full_list,]
pheatmap(heatmap_mat, scale = "row", color = pal, gaps_col = c(42, 64, 106, 112),
         cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T, annotation_col = sample_df, 
         annotation_names_col = T, annotation_legend = T, annotation_colors = anno_cols)

############################################################################################

# PT_only, PT_low, SDX_only, GS_low, GS_only, PT_XG_High, GS_SDX_High, PT_SDX_High
write.table(PT_only, "~/Documents/ProjectWorkspace/DEAnalysis/Annotation/GeneList_Henan/PT_only.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(PT_low, "~/Documents/ProjectWorkspace/DEAnalysis/Annotation/GeneList_Henan/PT_low.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(SDX_only, "~/Documents/ProjectWorkspace/DEAnalysis/Annotation/GeneList_Henan/SDX_only.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(GS_low, "~/Documents/ProjectWorkspace/DEAnalysis/Annotation/GeneList_Henan/GS_low.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(GS_only, "~/Documents/ProjectWorkspace/DEAnalysis/Annotation/GeneList_Henan/GS_only.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(PT_XG_High, "~/Documents/ProjectWorkspace/DEAnalysis/Annotation/GeneList_Henan/PT_XG_High.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(GS_SDX_High, "~/Documents/ProjectWorkspace/DEAnalysis/Annotation/GeneList_Henan/GS_SDX_High.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(PT_SDX_High, "~/Documents/ProjectWorkspace/DEAnalysis/Annotation/GeneList_Henan/PT_SDX_High.txt", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)








################################################################################################
# Depatched
################################################################################################

# pt_xg_resIHW$padj[is.na(pt_xg_resIHW$padj)] <- 1
# pt_gs_resIHW$padj[is.na(pt_gs_resIHW$padj)] <- 1
# xg_gs_resIHW$padj[is.na(xg_gs_resIHW$padj)] <- 1
# # filter for adjusted p less than 0.05 and lfc > 1
# pt_xg_sig <- pt_xg_resIHW[pt_xg_resIHW$padj < 0.05 & abs(pt_xg_resIHW$log2FoldChange) > 1,]
# pt_gs_sig <- pt_gs_resIHW[pt_gs_resIHW$padj < 0.05 & abs(pt_gs_resIHW$log2FoldChange) > 1,]
# xg_gs_sig <- xg_gs_resIHW[xg_gs_resIHW$padj < 0.05 & abs(xg_gs_resIHW$log2FoldChange) > 1,]
# # code genes by their pattern of differential expression
# all_sig_genes <- unique(c(rownames(pt_xg_sig), rownames(pt_gs_sig), rownames(xg_gs_sig)))
# check.gene <- function(gene, de_res){
#   if(gene %in% rownames(de_res)){
#     return(sign(de_res[rownames(de_res) == gene,]$stat))
#   } else {
#     return(0)
#   }
# }
# code_df <- data.frame(matrix(nrow = length(all_sig_genes), ncol = 3))
# rownames(code_df) <- all_sig_genes
# for(i in 1:length(all_sig_genes)){
#   temp_gene <- all_sig_genes[i]
#   a <- check.gene(temp_gene, pt_xg_sig)
#   b <- check.gene(temp_gene, pt_gs_sig)
#   c <- check.gene(temp_gene, xg_gs_sig)
#   code_df[i,] <- c(a,b,c)
# }
# # toString the codes
# string_code <- paste0(code_df[,1], code_df[,2], code_df[,3])
# names(string_code) <- all_sig_genes
# # write genes following DE patterns to separate files
# for(i in 1:length(unique(string_code))){
#   temp_code <- unique(string_code)[i]
#   temp_genes <- names(string_code[string_code == temp_code])
#   write.table(temp_genes, paste0(temp_code, "_DE_gene_list.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
# }
# 
# #
# setwd("~/Documents/ProjectWorkspace/DEAnalysis/Annotation/Reference/")
# 
# files <- list.files()
# GO_110_files <- files[grep("110_GO", files)]
# GO_110_genes <- c()
# for(i in 1:length(GO_110_files)){
#   temp <- read.table(GO_110_files[i], header = F, stringsAsFactors = F, sep = "\t")
#   print(length(setdiff(temp[,2], GO_110_genes)))
#   GO_110_genes <- unique(c(GO_110_genes, temp[,2]))
# }
# 
# GO_011_files <- files[grep("011_GO", files)]
# GO_011_genes <- c()
# for(i in 1:length(GO_011_files)){
#   temp <- read.table(GO_011_files[i], header = F, stringsAsFactors = F, sep = "\t")
#   print(length(setdiff(temp[,2], GO_011_genes)))
#   GO_011_genes <- unique(c(GO_011_genes, temp[,2]))
# }
# 
# # uaing x instead of dash because R doesn't allow dash in variable name
# GO_0x1x1_files <- files[grep("0-1-1_GO", files)]
# GO_0x1x1_genes <- c()
# for(i in 1:length(GO_0x1x1_files)){
#   temp <- read.table(GO_0x1x1_files[i], header = F, stringsAsFactors = F, sep = "\t")
#   print(length(setdiff(temp[,2], GO_0x1x1_genes)))
#   GO_0x1x1_genes <- unique(c(GO_0x1x1_genes, temp[,2]))
# }
# 
# # some GO genes have there names changed by GO, need to filter out for heatmap or map gene names
# filt_GO_110_genes <- GO_110_genes[GO_110_genes %in% rownames(assay(vsd))]
# filt_GO_011_genes <- GO_011_genes[GO_011_genes %in% rownames(assay(vsd))]
# filt_GO_0x1x1_genes <- GO_0x1x1_genes[GO_0x1x1_genes %in% rownames(assay(vsd))]
# 
# all_GO_genes <- c(filt_GO_110_genes[1:200], filt_GO_011_genes[1:150], filt_GO_0x1x1_genes[1:100])
# all_GO_genes_test <- c(filt_GO_110_genes, filt_GO_011_genes, filt_GO_0x1x1_genes)
# #m_ind <- match(all_GO_genes, rownames(filt_GO_110_genes))
# 
# #est <- estimateSizeFactors(new_dds)
# #mat <- log2((counts(est, normalized=TRUE)) + 1)
# #mat <- fpm(est)
# 
# heatmap_mat <- assay(vsd)[all_GO_genes,]
# #
# 
# #
# pal <- colorRampPalette(rev(brewer.pal("RdBu", n = 11)[c(1,2,3,4,6,8,9,10,11)]))(100)
# sample_df <- data.frame(Type = metadata$Type)
# rownames(sample_df) <- metadata$Sample
# gene_df <- data.frame(GO = c(rep("A", 200), rep("B", 150), rep("C", 100)))
# rownames(gene_df) <- all_GO_genes
# #
# pal_a <- brewer.pal("Dark2", n = 5)
# pal_b <- brewer.pal("Set1", n = 8)
# anno_cols <- list(Type = c("PT" = pal_a[1], "GS" = pal_a[2], "XG" = pal_a[3], "SDX" = pal_a[4], "XDS" = pal_a[5]), 
#                   GO   = c("A"  = pal_b[2], "B"  = pal_b[3], "C"  = pal_b[1]))
# #
# full_list <- c(PT_only, PT_SDX_High, PT_XG_High, GS_low, GS_only, PT_low, GS_SDX_High, SDX_only)
# back_list <- subset(rownames(assay(vsd)), !(rownames(assay(vsd)) %in% full_list))
# write.table(full_list, "~/Documents/ProjectWorkspace/DEAnalysis/Annotation/GeneList_Henan/full_list.txt", 
#             col.names = FALSE, row.names = FALSE, quote = FALSE)
# write.table(back_list, "~/Documents/ProjectWorkspace/DEAnalysis/Annotation/GeneList_Henan/back_list.txt", 
#             col.names = FALSE, row.names = FALSE, quote = FALSE)
# write.table(all_GO_genes, "~/Documents/ProjectWorkspace/DEAnalysis/Annotation/GeneList_Henan/all_GO_genes.Nick.txt", 
#             col.names = FALSE, row.names = FALSE, quote = FALSE)
# #
# pal <- colorRampPalette(rev(brewer.pal("RdBu", n = 11)[c(1,2,3,4,6,8,9,10,11)]))(100)
# sample_df <- data.frame(Type = metadata$Type)
# rownames(sample_df) <- metadata$Sample
# # PT_only, PT_low, SDX_only, GS_low, GS_only, PT_XG_High, GS_SDX_High, PT_SDX_High
# gene_df <- data.frame(Group = c(rep("PT_only"    , length(PT_only)), 
#                                 rep("PT_SDX_High", length(PT_SDX_High)),
#                                 rep("PT_XG_High" , length(PT_XG_High)),
#                                 rep("GS_low"     , length(GS_low)),
#                                 rep("GS_only"    , length(GS_only)),
#                                 rep("PT_low"     , length(PT_low)),
#                                 rep("GS_SDX_High", length(GS_SDX_High)),
#                                 rep("SDX_only"   , length(SDX_only))))
# rownames(gene_df) <- full_list
# #
# pal_a <- brewer.pal("Dark2", n = 5)
# pal_b <- brewer.pal("Set1", n = 8)
# anno_cols <- list(Type = c("PT" = pal_a[1], "GS" = pal_a[2], "XG" = pal_a[3], "SDX" = pal_a[4], "XDS" = pal_a[5]), 
#                   Group = c("PT_only"     = pal_b[1], 
#                             "PT_low"      = pal_b[2], 
#                             "SDX_only"    = pal_b[3],
#                             "GS_low"      = pal_b[4],
#                             "GS_only"     = pal_b[5],
#                             "PT_XG_High"  = pal_b[6],
#                             "GS_SDX_High" = pal_b[7],
#                             "PT_SDX_High" = pal_b[8]))
# #
# pheatmap(heatmap_mat, scale = "row", color = pal, gaps_col = c(42, 64, 106, 112),
#          cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T, 
#          annotation_row = gene_df, annotation_col = sample_df, 
#          annotation_names_row = T, annotation_names_col = T, 
#          annotation_legend = T, annotation_colors = anno_cols)
# 
# write.table(test_list, "~/Documents/ProjectWorkspace/DEAnalysis/Annotation/GeneList_Henan/test_list.txt", 
#             col.names = FALSE, row.names = FALSE, quote = FALSE)
# 
# #
# for (gene in all_genes) {
#   inter <- pt_gs_MA$data[grep(gene, pt_gs_MA$data$name),]$sig
#   inter <- pt_xg_MA$data[grep(gene, pt_xg_MA$data$name),]$sig
#   inter <- pt_sdx_MA$data[grep(gene, pt_sdx_MA$data$name),]$sig
#   inter <- xg_gs_MA$data[grep(gene, xg_gs_MA$data$name),]$sig
#   inter <- xg_sdx_MA$data[grep(gene, xg_sdx_MA$data$name),]$sig
#   inter <- gs_sdx_MA$data[grep(gene, gs_sdx_MA$data$name),]$sig
# }
# 
# #
# heatmap_mat <- assay(vsd)[all_GO_genes_test,]
# pheatmap(heatmap_mat, scale = "row", color = pal, gaps_col = c(42, 64, 106, 112),
#          cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T, annotation_col = sample_df, 
#          annotation_names_col = T, annotation_legend = T, annotation_colors = anno_cols)
# 
# 
# filt_GO_110_genes[1:200], filt_GO_011_genes[1:150], filt_GO_0x1x1_genes[1:100]
#
# gene_list <- lapply(1:length(all_GO_genes), function(x) matrix(rep(0, 16), nrow=4, ncol=4, dimnames = list(c("PT", "XG", "GS", "SDX"), c("PT", "XG", "GS", "SDX"))))
# names(gene_list) <- all_GO_genes
# ################################ fail to identify #######################################
# matrix <- rbind(c(0,0,0,0),
#                 c(0,0,0,0),
#                 c(0,0,0,0),
#                 c(0,0,0,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# none_only <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     none_only <- c(none_only, gene)
#   }
# }
# 
# ################################ Only high expressd in PT ################################
# matrix <- rbind(c(0,1,1,1),
#                 c(0,0,0,0),
#                 c(0,0,0,0),
#                 c(0,0,0,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# PT_only <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     PT_only <- c(PT_only, gene)
#   }
# }
# 
# ################################  Only high expressd in XG ################################ 
# matrix <- rbind(c(0,0,0,0),
#                 c(1,0,1,1),
#                 c(0,0,0,0),
#                 c(0,0,0,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# # Only high expressd in PT
# xg_only <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     xg_only <- c(xg_only, gene)
#   }
# }
# 
# ################################  Only high expressd in GS ################################ 
# matrix <- rbind(c(0,0,0,0),
#                 c(0,0,0,0),
#                 c(1,1,0,1),
#                 c(0,0,0,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# GS_only <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     GS_only <- c(GS_only, gene)
#   }
# }
# 
# ################################ Only high expressd in SDX ################################ 
# matrix <- rbind(c(0,0,0,0),
#                 c(0,0,0,0),
#                 c(0,0,0,0),
#                 c(1,1,1,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# SDX_only <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     SDX_only <- c(SDX_only, gene)
#   }
# }
# 
# ################################ Only low expressd in PT ################################ 
# matrix <- rbind(c(0,0,0,0),
#                 c(1,0,0,0),
#                 c(1,0,0,0),
#                 c(1,0,0,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# PT_low <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     PT_low <- c(PT_low, gene)
#   }
# }
# 
# ################################ Only low expressd in XG ################################ 
# matrix <- rbind(c(0,1,0,0),
#                 c(0,0,0,0),
#                 c(0,1,0,0),
#                 c(0,1,0,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# XG_low <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     XG_low <- c(XG_low, gene)
#   }
# }
# 
# ################################ Only low expressd in GS ################################ 
# matrix <- rbind(c(0,0,1,0),
#                 c(0,0,1,0),
#                 c(0,0,0,0),
#                 c(0,0,1,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# GS_low <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     GS_low <- c(GS_low, gene)
#   }
# }
# 
# ################################ Only low expressd in SDX ################################ 
# matrix <- rbind(c(0,0,0,1),
#                 c(0,0,0,1),
#                 c(0,0,0,1),
#                 c(0,0,0,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# SDX_low <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     SDX_low <- c(SDX_low, gene)
#   }
# }
# 
# ################################ High express in PT and XG ################################ 
# matrix <- rbind(c(0,0,1,1),
#                 c(0,0,1,1),
#                 c(0,0,0,0),
#                 c(0,0,0,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# PT_XG_High <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     PT_XG_High <- c(PT_XG_High, gene)
#   }
# }
# 
# ################################ High express in PT and GS ################################ 
# matrix <- rbind(c(0,1,0,1),
#                 c(0,0,0,0),
#                 c(0,1,0,1),
#                 c(0,0,0,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# PT_GS_High <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     PT_GS_High <- c(PT_GS_High, gene)
#   }
# }
# 
# ################################ High express in PT and SDX ################################ 
# matrix <- rbind(c(0,1,1,0),
#                 c(0,0,0,0),
#                 c(0,0,0,0),
#                 c(0,1,1,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# PT_SDX_High <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     PT_SDX_High <- c(PT_SDX_High, gene)
#   }
# }
# 
# ################################ High express in XG and GS ################################ 
# matrix <- rbind(c(0,0,0,0),
#                 c(1,0,0,1),
#                 c(1,0,0,1),
#                 c(0,0,0,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# XG_GS_High <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     XG_GS_High <- c(XG_GS_High, gene)
#   }
# }
# 
# ################################ High express in XG and SDX ################################ 
# matrix <- rbind(c(0,0,0,0),
#                 c(1,0,1,0),
#                 c(0,0,0,0),
#                 c(1,0,1,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# XG_SDX_High <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     XG_SDX_High <- c(XG_SDX_High, gene)
#   }
# }
# 
# ################################ High express in GS and SDX ################################ 
# matrix <- rbind(c(0,0,0,0),
#                 c(0,0,0,0),
#                 c(1,1,0,0),
#                 c(1,1,0,0))
# rownames(matrix) <- c("PT", "XG", "GS", "SDX")
# colnames(matrix) <- c("PT", "XG", "GS", "SDX")
# #
# GS_SDX_High <- c()
# for (gene in all_genes) {
#   if (identical(gene_list[[gene]], matrix)) {
#     GS_SDX_High <- c(GS_SDX_High, gene)
#   }
# }

# 
# pt_xg_resIHW   = readRDS("RDS/pt_xg_resIHW.rds"  )
# pt_gs_resIHW   = readRDS("RDS/pt_gs_resIHW.rds"  )
# pt_sdx_resIHW  = readRDS("RDS/pt_sdx_resIHW.rds" )
# pt_xds_resIHW  = readRDS("RDS/pt_xds_resIHW.rds" )
# xg_gs_resIHW   = readRDS("RDS/xg_gs_resIHW.rds"  )
# xg_sdx_resIHW  = readRDS("RDS/xg_sdx_resIHW.rds" )
# xg_xds_resIHW  = readRDS("RDS/xg_xds_resIHW.rds" )
# gs_sdx_resIHW  = readRDS("RDS/gs_sdx_resIHW.rds" )
# gs_xds_resIHW  = readRDS("RDS/gs_xds_resIHW.rds" )
# sdx_xds_resIHW = readRDS("RDS/sdx_xds_resIHW.rds")

#
# for (gene in all_GO_genes) {
#   # 
#   inter <- pt_gs_MA$data[grep(gene, pt_gs_MA$data$name),]$sig
#   if (length(inter) != 0) {
#     if (grepl(pattern = "Up", x = inter)) {
#       gene_list[[gene]]["PT", "GS"] <- 1
#     } else if (grepl(pattern = "Down", x = inter)) {
#       gene_list[[gene]]["GS", "PT"] <- 0
#     }
#   }
#   #
#   inter <- pt_xg_MA$data[grep(gene, pt_xg_MA$data$name),]$sig
#   if (length(inter) != 0) {
#     if (grepl(pattern = "Up", x = inter)) {
#       gene_list[[gene]]["PT", "XG"] <- 1
#     } else if (grepl(pattern = "Down", x = inter)) {
#       gene_list[[gene]]["XG", "PT"] <- 0
#     }
#   }
#   #
#   inter <- pt_sdx_MA$data[grep(gene, pt_sdx_MA$data$name),]$sig
#   if (length(inter) != 0) {
#     if (grepl(pattern = "Up", x = inter)) {
#       gene_list[[gene]]["PT", "SDX"] <- 1
#     } else if (grepl(pattern = "Down", x = inter)) {
#       gene_list[[gene]]["SDX", "PT"] <- 0
#     }
#   }
#   #
#   inter <- xg_gs_MA$data[grep(gene, xg_gs_MA$data$name),]$sig
#   if (length(inter) != 0) {
#     if (grepl(pattern = "Up", x = inter)) {
#       gene_list[[gene]]["XG", "GS"] <- 1
#     } else if (grepl(pattern = "Down", x = inter)) {
#       gene_list[[gene]]["GS", "XG"] <- 0
#     }
#   }
#   #
#   inter <- xg_sdx_MA$data[grep(gene, xg_sdx_MA$data$name),]$sig
#   if (length(inter) != 0) {
#     if (grepl(pattern = "Up", x = inter)) {
#       gene_list[[gene]]["XG", "SDX"] <- 1
#     } else if (grepl(pattern = "Down", x = inter)) {
#       gene_list[[gene]]["SDX", "XG"] <- 0
#     }
#   }
#   #
#   inter <- gs_sdx_MA$data[grep(gene, gs_sdx_MA$data$name),]$sig
#   if (length(inter) != 0) {
#     if (grepl(pattern = "Up", x = inter)) {
#       gene_list[[gene]]["GS", "SDX"] <- 1
#     } else if (grepl(pattern = "Down", x = inter)) {
#       gene_list[[gene]]["SDX", "GS"] <- 0
#     }
#   }
# }