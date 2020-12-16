#############################################################################################################################################
#
#   WIREs_RNA_Review.R
#
#   Perform scRNA-seq analysis on cells from Loo et al., 2019 (Nature Communications) for WIREs RNA Review.
#
#   
#
#   Nicholas Page, December 2020
#   Rasin Lab, Dept. Neuroscience and Cell Biology, Rutgers-New Brunswick
#   Neuroscience Graduate Program, University of California, San Francisco
#
#############################################################################################################################################

rm(list=ls())
options(stringsAsFactors = FALSE)

library(gridExtra)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(Seurat)
library(Cairo)
library(dplyr)
library(EWCE)

setwd('C:/Users/npage/Desktop/WIREs RNA Review/')

#############################################################################################################################################
#
#   Pipeline:
#
#   (1) Download scRNA-seq data from Loo et al and perform gene and cell filtering
#   (2) Recreate UMAP plots for cells retained in this analysis based on original clusters from Loo et al
#   (3) Filter to cell-subtypes of interest and then creat UMAP plots fo meta cell-types
#   (4) Create cell-type specificity matrices with EWCE
#   (5) Analyze specificity results and determine top ranked RBPs
#   (6) Plot genes of interest and prepare figure
#
#############################################################################################################################################

### (1) Download scRNA-seq data from Loo et al and perform gene and cell filtering ##########################################################

# load E14.5.5 cells and filter genes that have a total expression of less than 50 across all cells and 
# filter cells that have less than or equal to 1000 total features (genes) detected
raw_data = read.table('./GSE123335_E14_combined_matrix.txt', row.names = 1, header = TRUE, sep = '\t')
E14.5 <- CreateSeuratObject(counts = raw_data, min.cells = 20, min.features = 1000, project = "E14.5")

E14.5_cell_types = read.table('./GSE123335_E14_combined_matrix_ClusterAnnotations.txt', row.names = 1, header = TRUE, sep = '\t')

# Filter E14.5 dataset to cells with total counts within half a standard deviation of the mean
E14.5 <- subset(E14.5, subset = nCount_RNA >= mean(E14.5$nCount_RNA) - 0.5*sd(E14.5$nCount_RNA) & 
                  nCount_RNA <= mean(E14.5$nCount_RNA) + 0.5*sd(E14.5$nCount_RNA))

FeatureScatter(object = E14.5, feature1 = "nFeature_RNA", feature2 = "nCount_RNA") ### Looks good!!



# load P0 cells and filter genes that have a total expression of less than 50 across all cells and 
# filter cells that have less than or equal to 1000 total features (genes) detected
raw_data = read.table('./GSE123335_P0_combined_matrix.txt', row.names = 1, header = TRUE, sep = '\t')
P0 <- CreateSeuratObject(counts = raw_data, min.cells = 20, min.features = 1000, project = "P0")

P0_cell_types = read.table('./GSE123335_P0_combined_matrix_ClusterAnnotations.txt', row.names = 1, header = TRUE, sep = '\t')

# Filter P0 dataset to cells with total counts within half a standard deviation of the mean
P0 <- subset(P0, subset = nCount_RNA >= mean(P0$nCount_RNA) - 0.5*sd(P0$nCount_RNA) & 
                  nCount_RNA <= mean(P0$nCount_RNA) + 0.5*sd(P0$nCount_RNA))

FeatureScatter(object = P0, feature1 = "nFeature_RNA", feature2 = "nCount_RNA") ### Looks good!!



### (2) Recreate UMAP plots for cells retained in this analysis based on original clusters from Loo et al ###################################

# log normalize the expression data and multiply by a scaling factor of 10000
E14.5 <- NormalizeData(object = E14.5, normalization.method = "LogNormalize", scale.factor = 10000)

# find genes that represent 'highly variable features' that exhibit high cell-to-cell variation within the dataset using the 'vst'
# method (see '?FindVariableFeatures') and then run PCA on the top 2000 variable features detected
E14.5 <- FindVariableFeatures(object = E14.5, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x = E14.5)
E14.5 <- ScaleData(object = E14.5, features = all.genes)
E14.5 <- RunPCA(object = E14.5, features = VariableFeatures(object = E14.5))

# perform t-SNE on E14.5 samples
E14.5 <- RunUMAP(object = E14.5, dims = 1:30)

# assign cell-type identities from Loo et al to cells in this analysis
E14.5_cell_types <- data.frame("Cell" = rownames(E14.5_cell_types), "Cluster" = E14.5_cell_types$Cluster)
E14.5_cell_types <- E14.5_cell_types[which(E14.5_cell_types$Cell %in% names(E14.5$orig.ident)),]
Idents(object = E14.5) <- E14.5_cell_types$Cluster[match(E14.5_cell_types$Cell, names(E14.5$orig.ident))]
E14.5$orig.ident <- Idents(object = E14.5)

DimPlot(object = E14.5, reduction = "umap", pt.size = 2) + ggtitle('E14.5')



# log normalize the expression data and multiply by a scaling factor of 10000
P0 <- NormalizeData(object = P0, normalization.method = "LogNormalize", scale.factor = 10000)

# find genes that represent 'highly variable features' that exhibit high cell-to-cell variation within the dataset using the 'vst'
# method (see '?FindVariableFeatures') and then run PCA on the top 2000 variable features detected
P0 <- FindVariableFeatures(object = P0, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x = P0)
P0 <- ScaleData(object = P0, features = all.genes)
P0 <- RunPCA(object = P0, features = VariableFeatures(object = P0))

# perform t-SNE on P0 samples
P0 <- RunUMAP(object = P0, dims = 1:30)

# assign cell-type identities from Loo et al to cells in this analysis
P0_cell_types <- data.frame("Cell" = rownames(P0_cell_types), "Cluster" = P0_cell_types$Cluster)
P0_cell_types <- P0_cell_types[which(P0_cell_types$Cell %in% names(P0$orig.ident)),]
P0 <- subset(P0, cells = WhichCells(P0, which((rownames(P0@meta.data) %in% P0_cell_types$Cell) == TRUE)))
Idents(object = P0) <- P0_cell_types$Cluster[match(P0_cell_types$Cell, names(P0$orig.ident))]
P0$orig.ident <- Idents(object = P0)

DimPlot(object = P0, reduction = "umap", pt.size = 2) + ggtitle('P0')



### (3) Filter to cell-subtypes of interest and then creat UMAP plots fo meta cell-types ####################################################

# filter cell-subtypes of interest for the analysis
E14.5_subtypes <- subset(E14.5, idents = c("RG1 [8-E]",
                                     "LayerV-VI [3-E]",
                                     "LayerV-VI [5-E]",
                                     "Int2 [12-E]",
                                     "LayerV-VI [7-E]",
                                     "RG4 [10-E]",
                                     "SVZ3 (proliferating) [15-E]",
                                     "SVZ2 (VZ-SVZ) [11-E]",
                                     "SVZ1 (migrating) [4-E]" ,
                                     "LayerV-VI [13-E]",
                                     "Int1 [1-E]",
                                     "RG2 [14-E]",
                                     "LayerV-VI [2-E]"))

# Reassign cells to simplified cell-types
new.cluster.ids <- c('RG', 'ExN', 'ExN', 'InN', 'ExN', 'RG', 'IPC', 'IPC', 'IPC', 'ExN', 'InN', 'RG', 'ExN')
names(new.cluster.ids) <- levels(E14.5_subtypes)
E14.5_types <- RenameIdents(E14.5_subtypes, new.cluster.ids)

# Replot UMAP with simplified cell-types
DimPlot(object = E14.5_types, reduction = "umap", pt.size = 2, group.by = "ident") + ggtitle('E14.5')



# filter cell-subtypes of interest for the analysis
P0_subtypes <- subset(P0, idents = c("Int1 [5-P]",
                                   "LayerII-IV [1-P]",
                                   "LayerII-IV [4-P]",
                                   "Oligodendrocytes [16-P]",
                                   "Astrocytes (immature) 1 [10-P]",
                                   "Int4 [6-P]",
                                   "Astrocytes (immature) 2 [13-P]",
                                   "Int2 [14-P]" ,
                                   "Layer II-IV [15-P]",
                                   "LayerV-VI [12-P]",
                                   "Int3 [11-P]",
                                   "LayerV-VI [18-P]"))

# Reassign cells to simplified cell-types
new.cluster.ids <- c('InN', 'ExN', 'ExN', 'Oligo', 'Astro', 'InN', 'Astro', 'InN', 'ExN', 'ExN', 'InN', 'ExN')
names(new.cluster.ids) <- levels(P0_subtypes)
P0_types <- RenameIdents(P0_subtypes, new.cluster.ids)

# Replot UMAP with simplified cell-types
DimPlot(object = P0_types, reduction = "umap", pt.size = 2, group.by = "ident") + ggtitle('P0')



### (4) Create cell-type specificity matrices with EWCE #####################################################################################

E14.5_expr_matrix <- GetAssayData(object = E14.5_types, slot = "counts")
E14.5_cell_type_assignments <- cbind(names(Idents(E14.5_types)), as.vector(Idents(E14.5_types)), as.vector(E14.5_types@meta.data$orig.ident))

### Double checks that the order of the samples is the same as the metadata
### THIS STEP IS VERY IMPORTANT!!!
print(identical(E14.5_cell_type_assignments[,1], colnames(E14.5_expr_matrix)))

colnames(E14.5_cell_type_assignments) <- c('cell_id', 'level1class', 'level2class')
E14.5_cell_type_assignments <- data.frame(E14.5_cell_type_assignments)

### convert expr matrix ###
E14.5_expr_matrix <- as.matrix(E14.5_expr_matrix)

# Generate cell-type data for just the cortex
E14.5_exp_CortexOnly_DROPPED = drop.uninformative.genes(exp=E14.5_expr_matrix, level2annot = E14.5_cell_type_assignments$level2class)

# Generate specificity matric
E14.5_annotLevels = list(level1class=E14.5_cell_type_assignments$level1class,level2class=E14.5_cell_type_assignments$level2class)
E14.5_fNames_CortexOnly = generate.celltype.data(exp=E14.5_exp_CortexOnly_DROPPED,annotLevels=E14.5_annotLevels,groupName="E14.5_kiCortexOnly")
load(paste0('./CellTypeData_E14.5_kiCortexOnly.rda'))

# Check specificity matrix
E14.5_specificity <- data.frame(ctd[[1]]$specificity)
head(E14.5_specificity)



P0_expr_matrix <- GetAssayData(object = P0_types, slot = "counts")
P0_cell_type_assignments <- cbind(names(Idents(P0_types)), as.vector(Idents(P0_types)), as.vector(P0_types@meta.data$orig.ident))

### Double checks that the order of the samples is the same as the metadata
### THIS STEP IS VERY IMPORTANT!!!
print(identical(P0_cell_type_assignments[,1], colnames(P0_expr_matrix)))

colnames(P0_cell_type_assignments) <- c('cell_id', 'level1class', 'level2class')
P0_cell_type_assignments <- data.frame(P0_cell_type_assignments)

### convert expr matrix ###
P0_expr_matrix <- as.matrix(P0_expr_matrix)

# Generate cell-type data for just the cortex
P0_exp_CortexOnly_DROPPED = drop.uninformative.genes(exp=P0_expr_matrix, level2annot = P0_cell_type_assignments$level2class)

# Generate specificity matric
P0_annotLevels = list(level1class=P0_cell_type_assignments$level1class,level2class=P0_cell_type_assignments$level2class)
P0_fNames_CortexOnly = generate.celltype.data(exp=P0_exp_CortexOnly_DROPPED,annotLevels=P0_annotLevels,groupName="P0_kiCortexOnly")
load(paste0('./CellTypeData_P0_kiCortexOnly.rda'))

# Check specificity matrix
P0_specificity <- data.frame(ctd[[1]]$specificity)
head(P0_specificity)



### (5) Analyze specificity results and determine top ranked RBPs ###########################################################################

# Get the max specificity for each gene and order by this metrix
E14.5_specificity$max <- as.vector(apply(E14.5_specificity, 1, max))
E14.5_specificity <- E14.5_specificity[order(E14.5_specificity$max, decreasing = TRUE),]
head(E14.5_specificity, 100)

# Get the max specificity for each gene and order by this metrix
P0_specificity$max <- as.vector(apply(P0_specificity, 1, max))
P0_specificity <- P0_specificity[order(P0_specificity$max, decreasing = TRUE),]
head(P0_specificity, 100)


# get list of TFs and RBPs 
TFs <- read.table(file = './GO_term_summary_20200722_205508.csv', header = TRUE, sep = ',')
TFs <- unique(TFs$Symbol)

RBPs <- read.table(file = './GO_term_summary_20200722_205617.csv', header = TRUE, sep = ',')
RBPs <- unique(RBPs$Symbol)

# Set up for specificity analysis
E14.5_specificity$rank <- seq.int(nrow(E14.5_specificity))
E14.5_specificity$TFs <- ifelse(rownames(E14.5_specificity) %in% TFs, E14.5_specificity$rank, 0)
E14.5_specificity$RBPs <- ifelse(rownames(E14.5_specificity) %in% RBPs, E14.5_specificity$rank, 0)
E14.5_specificity$GO <- ifelse(E14.5_specificity$TFs > E14.5_specificity$RBPs, 'TF', 
                         ifelse(E14.5_specificity$RBPs > E14.5_specificity$TFs, 'RBP', 'NA'))
E14.5_specificity$Rank <- ifelse(E14.5_specificity$TFs > E14.5_specificity$RBPs, E14.5_specificity$TFs, 
                           ifelse(E14.5_specificity$RBPs > E14.5_specificity$TFs, E14.5_specificity$RBPs, 0))
E14.5_to_plot <- data.frame(E14.5_specificity$Rank, E14.5_specificity$GO)
rownames(E14.5_to_plot) <- rownames(E14.5_specificity)
colnames(E14.5_to_plot) <- c('Rank', 'GO')
E14.5_to_plot <- E14.5_to_plot[which(E14.5_to_plot$GO != 'NA'),]
head(E14.5_to_plot)

# Create plot
E14.5_density_plot <- ggdensity(E14.5_to_plot, x = "Rank",
                             add = "median", rug = TRUE,
                             color = "GO", fill = "GO",
                             palette = c("#0073C2FF", "#FC4E07")) +
  scale_x_continuous(breaks = c(0, 4000, 8000, 12000)) +
  ylab("Density") +
  theme(plot.title = element_blank(), 
        axis.title.x = element_text(size = 10, margin = margin(t = 2, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 2, b = 0, l = 0)),
        axis.text.y.left = element_blank(),
        legend.text=element_text(size = 8),
        legend.title = element_blank(),
        legend.position = 'right',
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.spacing.y = unit(1, "mm"),
        legend.box.margin = margin(-10, -5, -10, -10),
        plot.margin = margin(15, 5, 15, 15))
plot(E14.5_density_plot)

# Examine RBP List
E14.5_RBP_specificity <- E14.5_specificity[which(E14.5_specificity$GO == "RBP"),]
head(E14.5_RBP_specificity, 50)

# Examine TF List
E14.5_TF_specificity <- E14.5_specificity[which(E14.5_specificity$GO == "TF"),]
head(E14.5_TF_specificity, 50)



# Set up for specificity analysis
P0_specificity$rank <- seq.int(nrow(P0_specificity))
P0_specificity$TFs <- ifelse(rownames(P0_specificity) %in% TFs, P0_specificity$rank, 0)
P0_specificity$RBPs <- ifelse(rownames(P0_specificity) %in% RBPs, P0_specificity$rank, 0)
P0_specificity$GO <- ifelse(P0_specificity$TFs > P0_specificity$RBPs, 'TF', 
                               ifelse(P0_specificity$RBPs > P0_specificity$TFs, 'RBP', 'NA'))

P0_specificity$Rank <- ifelse(P0_specificity$TFs > P0_specificity$RBPs, P0_specificity$TFs, 
                                 ifelse(P0_specificity$RBPs > P0_specificity$TFs, P0_specificity$RBPs, 0))
P0_to_plot <- data.frame(P0_specificity$Rank, P0_specificity$GO)
rownames(P0_to_plot) <- rownames(P0_specificity)
colnames(P0_to_plot) <- c('Rank', 'GO')
P0_to_plot <- P0_to_plot[which(P0_to_plot$GO != 'NA'),]
head(P0_to_plot)

# Create plot
P0_density_plot <- ggdensity(P0_to_plot, x = "Rank",
          add = "median", rug = TRUE,
          color = "GO", fill = "GO",
          palette = c("#0073C2FF", "#FC4E07")) +
  scale_x_continuous(breaks = c(0, 4000, 8000, 12000)) +
  ylab("Density") +
  theme(plot.title = element_blank(), 
        axis.title.x = element_text(size = 10, margin = margin(t = 2, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 2, b = 0, l = 0)),
        axis.text.y.left = element_blank(),
        legend.text=element_text(size = 8),
        legend.title = element_blank(),
        legend.position = 'right',
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.spacing.y = unit(1, "mm"),
        legend.box.margin = margin(-10, -5, -10, -10),
        plot.margin = margin(15, 5, 15, 15))
plot(P0_density_plot)

# Examine RBP List
P0_RBP_specificity <- P0_specificity[which(P0_specificity$GO == "RBP"),]
head(P0_RBP_specificity, 50)

# Examine TF List
P0_TF_specificity <- P0_specificity[which(P0_specificity$GO == "TF"),]
head(P0_TF_specificity, 50)

### (6) Plot genes of interest and prepare figure ###########################################################################################

E14.5_UMAP <- DimPlot(object = E14.5_types, reduction = "umap", pt.size = 0.5) + 
  labs(color = 'E14.5') + theme_classic() + 
  theme(plot.title = element_blank(), 
                     axis.title.x = element_text(size = 10, margin = margin(t = 2, r = 0, b = 0, l = 0)),
                     axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 2, b = 0, l = 0)),
                     legend.text=element_text(size = 8),
                     legend.position = 'right',
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     strip.text.x = element_text(size = 10),
                     legend.spacing.y = unit(1, "mm"),
                     legend.box.margin = margin(-10, -5, -10, -10),
                     plot.margin = margin(15, 15, 15, 15))

P0_UMAP <- DimPlot(object = P0_types, reduction = "umap", pt.size = 0.5) + 
  labs(color = "P0") + theme_classic() + 
  xlim(-8.5, 8.5) + ylim(-10, 13) +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 10, margin = margin(t = 2, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 2, b = 0, l = 0)),
        legend.text=element_text(size = 8),
        legend.position = 'right',
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size = 10),
        legend.spacing.y = unit(1, "mm"),
        legend.box.margin = margin(-10, -5, -10, -10),
        plot.margin = margin(15, 15, 15, 15))

# Creat graph for specificity violin vln_plot
P0_to_vln_plot <- data.frame(P0_specificity$max, P0_specificity$GO)
rownames(P0_to_vln_plot) <- rownames(P0_specificity)
colnames(P0_to_vln_plot) <- c('Specificity', 'GO')
P0_to_vln_plot <- P0_to_vln_plot[which(P0_to_vln_plot$GO != 'NA'),]
head(P0_to_vln_plot)

# Create plot
P0_vln_plot <- ggplot(P0_to_vln_plot, aes(x = GO, y = Specificity, color = GO, fill = GO)) +
  geom_point(P0_to_vln_plot[which(rownames(P0_to_vln_plot) %in% c("Neurod2", "Dlx5", 'Sox10', "Gli2")),], 
             mapping = aes(x = GO, y = Specificity), inherit.aes = FALSE, size = 2) +
  geom_text_repel(P0_to_vln_plot[which(rownames(P0_to_vln_plot) %in% c("Neurod2", "Dlx5", 'Sox10', "Gli2")),], 
                  mapping = aes(x = GO, y = Specificity, label = rownames(P0_to_vln_plot[which(rownames(P0_to_vln_plot) %in% c("Neurod2", "Dlx5", 'Sox10', "Gli2")),])), inherit.aes = FALSE) +
  geom_violin() +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1)) + 
  scale_color_manual(values = c("#0073C2FF", "#FC4E07")) +
  scale_fill_manual(values = alpha(c("#0073C2FF", "#FC4E07"), 0.1)) +
  ylab("Specificity Score") +
  geom_violin(size = 1) +
  theme_classic() +
  theme(plot.title = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 2, b = 0, l = 0)),
        legend.text=element_text(size = 8),
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.spacing.y = unit(1, "mm"),
        legend.box.margin = margin(-10, -5, -10, -10),
        plot.margin = margin(15, 5, 15, 15))
plot(P0_vln_plot)


Zfp36l1_plot <- FeaturePlot(object = E14.5_types, reduction = 'umap', features = c('Zfp36l1'), pt.size = 0.2, cols = c('gray95', 'darkred')) + 
  ggtitle('E14.5 Zfp36l1') +
  theme_classic() +
  theme(axis.title.x = element_text(size = 10, margin = margin(t = 0, r = 2, b = 0, l = 0)),
        axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 2, b = 0, l = 0)),
        legend.text=element_text(size = 8),
        legend.title = element_blank(),
        legend.position = 'right',
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size = 10),
        legend.spacing.y = unit(0.5, "mm"),
        legend.box.margin = margin(-10, -5, -10, -10),
        plot.margin = margin(20, 20, 20, 20))
plot(Zfp36l1_plot)

Zfp36l2_plot <- FeaturePlot(object = E14.5_types, reduction = 'umap', features = c('Zfp36l2'), pt.size = 0.2, cols = c('gray95', 'darkred')) + 
  ggtitle('E14.5 Zfp36l2') +
  theme_classic() +
  theme(axis.title.x = element_text(size = 10, margin = margin(t = 0, r = 2, b = 0, l = 0)),
        axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 2, b = 0, l = 0)),
        legend.text=element_text(size = 8),
        legend.title = element_blank(),
        legend.position = 'right',
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_text(size = 10),
        legend.spacing.y = unit(0.5, "mm"),
        legend.box.margin = margin(-10, -5, -10, -10),
        plot.margin = margin(20, 20, 20, 20))
plot(Zfp36l2_plot)



g1 <- ggplot()
plot_list <- list(E14.5_UMAP,
                  g1,
                  E14.5_density_plot,
                  P0_UMAP,
                  g1,
                  P0_density_plot,
                  P0_vln_plot,
                  Zfp36l1_plot, 
                  Zfp36l2_plot)

layout_matrix <- rbind(c(1, 1, 2, 2, 3, 3),
                       c(4, 4, 5, 5, 6, 6),
                       c(7, 7, 8, 8, 9, 9))

cairo_pdf('./Fig 1 Draft.pdf', width = 7.5, height = 6)
grid.arrange(grobs = plot_list, 
             layout_matrix = layout_matrix,
             padding = unit(2, 'line'))
dev.off()





