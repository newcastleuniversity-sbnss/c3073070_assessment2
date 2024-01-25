### Case study 3 - practical 8 ###
### RNA sequencing ###

##Install required packages
install.packages('BiocManager')
library(BiocManager)
install(c('tximport', 'DESeq2', 'biomaRt', 'pheatmap'))
install.packages('RColorBrewer')
install.packages('ggrepel')
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("viridisLite")
install.packages("dendextend")

##Library

library(tximport)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(tibble)
library(ggrepel)
library(viridisLite)
library(dendextend)


### 1.0 Download data ###
sample_table = read_csv('https://raw.githubusercontent.com/sjcockell/mmb8052/main/practicals/practical_08/data/sample_table.csv')
files = pull(sample_table, Run)
files = paste0('counts/', files, '/quant.sf')
names(files) = pull(sample_table, Run)
gene_map = read_csv('https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/extdata/gene_map.csv')
txi = tximport(files,
               type='salmon',
               tx2gene=gene_map,
               ignoreTxVersion=TRUE)

class(gene_map)



###  2.0 Differential gene expression analysis using DESeq2 ###

dds = DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ Group)
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds)
# the above 3 commands can be replaced by dds = DESeq(dds)
dds = DESeq(dds)
plotDispEsts(dds)


disp_plot <- plotDispEsts(dds, main = "Dispersion Estimates",
                          genecol = "black",
                          fitcol = "red",
                          finalcol = "cornflowerblue",
                          xlab = "Mean of nomalised counts",
                          ylab = "Dispersion")

disp_plot





### 3.0 Data Control ###
## transform data to log2 scale ##
#minimises differences between samples with small counts
#normalises with respect to library size

rld = rlog(dds)

##Scree plot
#determine the number of principal components to use for pca plot

sample_distance = dist(t(assay(rld)), method='euclidian')
sample_distance_matrix = as.matrix(sample_distance)


# eigenvalues determine the amount of variation explained by each principal component
pca_result2 <- prcomp(sample_distance_matrix, scale = TRUE)
eigenvalues <- pca_result2$sdev^2
line_colours <- viridis(length(eigenvalues), option = "plasma")


plot(1:length(eigenvalues), eigenvalues, type = "b", pch = 16,
     xlab = "Dimensions",
     ylab = "Eigenvalue",
     main = "Scree Plot",
     col = line_colours)

abline(h = 1, col = "indianred2", lty = 2)
axis(1, at = 1:length(eigenvalues), labels = seq(1, length(eigenvalues)))

text(1:length(eigenvalues), eigenvalues, labels = round(eigenvalues, 2), pos = 1, col = "black", cex = 0.8)

# Conduct principal component analysis
# Initial PCA plot
pca_data <- plotPCA(rld, intgroup='Group', returnData = TRUE)
pca_data


#PCA ggplot version to customsise pca
pca_1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA for alveolar macrophages under three treatment conditions",
       x = "PC1: 60% variance",
       y = "PC2: 18% variance") +
  scale_color_manual(values = c("Naive" = "#E69F00",
                                "Allo2h" = "#CC79A7",
                                "Allo24h" = "#56B4E9")) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

#ellipses added
pca_1 <- pca_1 + stat_ellipse(geom = "polygon", aes(group = Group, color = Group, fill = Group), alpha = 0.2) +
  scale_color_manual(values = c("Naive" = "#E69F00", "Allo2h" = "#CC79A7", "Allo24h" = "#56B4E9")) +
  scale_shape_manual(values = c("Naive" = 19, "Allo2h" = 17, "Allo24h" = 15)) +
  scale_fill_manual(values = c("Naive" = "#E69F00", "Allo2h" = "#CC79A7", "Allo24h" = "#56B4E9"))
pca_1

## 3d pca plot for three dimnensional visualisation ##
## NOTE: 3D PCA was not utilised for this dataet as the eigenvalues determined 2 PCs to explain enough variance for a 2D PCA. 
        #This script can be uncommented for use in datasets in which 3 principal components are necessary  to cover sufficient variance.

#pca_result <- prcomp(t(assay(rld)))

# Extract PC scores
#pca_scores <- as.data.frame(pca_result$x)

#print(pca_scores)

# Include Group information
#pca_scores$Group <- colData(rld)$Group
#pca_scores$Sample_Name <- colData(rld)$Sample_Name
#str(rld)
#view(pca_scores)
#colnames(colData(rld))

#pca_scores$Group <- factor(pca_scores$Group, levels = c("Naive", "Allo2h", "Allo24h"))
#color_scale <- c("Naive" = "#E69F00", "Allo2h" = "#CC79A7", "Allo24h" = "#56B4E9")

#pca_plot_3d <- plot_ly(
#data = pca_scores,
#x = ~PC1,
#y = ~PC2,
#z = ~PC3,
#color = ~Group,
#colors= color_scale,
#shape = ~Group,
# text = ~Sample_Name,
#type = "scatter3d",
#mode = "markers")


#layout <- list(
#title = "3D PCA Plot",
#scene = list(
#xaxis = list(title = "PC1"),
#yaxis = list(title = "PC2"),
#zaxis = list(title = "PC3")))

#pca_plot_3d <- pca_plot_3d %>% layout(layout)
#pca_plot_3d
#str(pca_scores)




## heirarchical clustering ##
# Dendrogram

heatmap_annotation <- data.frame(group = colData(dds)[, c('Group')], row.names = rownames(colData(dds)))

h <- hclust(sample_distance, method = "complete")
dend <- as.dendrogram(h)

dend %>%
  set("nodes_pch", 19) %>%
  set("nodes_cex", 0.7) %>%
  set("nodes_col", "black") %>%
  set("labels_col", value = c("Allo24h" = "#56B4E9", "Allo2h" = "#CC79A7", "Naive" = "#E69F00"), k=3) %>%
  set("branches_k_color", value = c("Allo24h" = "#56B4E9", "Allo2h" = "#CC79A7", "Naive" = "#E69F00"), k=3) %>%
  plot(las = 1, main = "Hierarchical Clustering of samples")



legend("topright", legend = levels(heatmap_annotation$group),
       col = c("#56B4E9", "#CC79A7", "#E69F00"), pch = 19,
       title = "Sample Group")


order_dend <- order.dendrogram(dend)


# Ensure the correct column name is used in colData(dds)
heatmap_annotation = data.frame(group=colData(dds)[,c('Group')], row.names=rownames(colData(dds)))

print(heatmap_annotation)



annotation_colours <- list(
  group = c("Allo24h" = "#56B4E9",
            "Allo2h" = "#CC79A7",
            "Naive" = "#E69F00"))

#Heatmap

pheatmap(sample_distance_matrix,
         clustering_distance_rows=sample_distance,
         clustering_distance_cols=sample_distance,
         treeheight_row = 45,
         treeheight_col = 45,
         border_color = "white",
         fontsize_row = 8,
         fontsize_col = 8,
         main = "Sample distances for GSE116583",
         annotation_col = heatmap_annotation,
         annotation_row = heatmap_annotation,
         color = colorRampPalette(c("#009E73", "white"))(100),
         angle_col = 90,
         gaps_col = cumsum(c(5, 5)),
         annotation_colors = annotation_colours)

heatmap_annotation

print(heatmap_annotation)
print(sample_distance_matrix)


### 4.0 Exploration of differentially expressed genes ###
#filtering results to add a column for significance based on p value <0.05

#allo24 vs naive
results_table = results(dds, contrast= c('Group', 'Allo24h', 'Naive'))
summary(results_table)
# filter any rows with an NA anywhere:
results_tibble = as_tibble(results_table, rownames='ensembl_gene_id')
filtered_results = filter(results_tibble, complete.cases(results_tibble))
#Volcano plot (basic)
filtered_results = mutate(filtered_results, logPVal = -log10(padj))
filtered_results = mutate(filtered_results, significant=padj<0.05)
filtered_results

#allo2 vs naive ( same as above for this contrast)
results_table2 = results(dds, contrast= c('Group', 'Allo2h', 'Naive'))
summary(results_table2)

results_tibble2 = as_tibble(results_table2, rownames='ensembl_gene_id')
filtered_results2 = filter(results_tibble2, complete.cases(results_tibble2))

filtered_results2 = mutate(filtered_results2, logPVal = -log10(padj))
filtered_results2 = mutate(filtered_results2, significant=padj<0.05)


#add a column called combined_significance for identifying differentially expressed genes

filtered_results$combined_significance <- ifelse(
  filtered_results$significant == FALSE, 'non-significant',
  ifelse(filtered_results$log2FoldChange > 1, 'upregulated',
         ifelse(filtered_results$log2FoldChange < -1, 'downregulated', 'p value only')))

#assigning colours to category of significance
filtered_results$combined_significance <- factor(
  filtered_results$combined_significance,
  levels = c('downregulated', 'non-significant', 'upregulated', 'p value only'))


## Repeating above steps for 2h treatment group 
filtered_results2 = mutate(filtered_results2, significant=padj<0.05)
filtered_results2$combined_significance <- ifelse(
  filtered_results2$significant == FALSE, 'non-significant',
  ifelse(filtered_results2$log2FoldChange > 1, 'upregulated',
         ifelse(filtered_results2$log2FoldChange < -1, 'downregulated', 'p value only')))


filtered_results2$combined_significance <- factor(
  filtered_results2$combined_significance,
  levels = c('downregulated', 'non-significant', 'upregulated', 'p value only'))



## Specifying upper and lower axes limits for  logFC in volcano plots
# 24h treatment group

#lower
filtered_results %>%
  dplyr::select(log2FoldChange) %>%
  min() %>%
  floor()
# [1] -20

#upper
filtered_results %>%
  dplyr::select(log2FoldChange) %>%
  max() %>%
  ceiling()
# [1] 23

c(-20, 23) %>%
  abs() %>%
  max()
# [1] 23



# 2h treatment group
#lower 
filtered_results2 %>%
  dplyr::select(log2FoldChange) %>%
  min() %>%
  floor()
# [1] -9

#upper
filtered_results2 %>%
  dplyr::select(log2FoldChange) %>%
  max() %>%
  ceiling()
# [1] 24

c(-9, 24) %>%
  abs() %>%
  max()
# [1] 24


## Annotations
##Extracting Gnee annotations for the top 10 most differentially expressed genes using BioMart for volcanoplot labels

# 24h treatment group

ensembl108 = biomaRt::useEnsembl(biomart="ensembl", version=108)
ensembl108 = useDataset("mmusculus_gene_ensembl", mart=ensembl108)
annotation = getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                'start_position', 'end_position',
                                'strand', 'gene_biotype', 'external_gene_name',
                                'description'),
                   filters = 'ensembl_gene_id', values = filtered_results$ensembl_gene_id,
                   mart = ensembl108)
annot_results = left_join(filtered_results, annotation)
annot_results = arrange(annot_results, padj)
#View(head(annot_results, 10))
degs = filter(annot_results, abs(log2FoldChange) > 1 & padj < 0.05)
top_10_DE <- (head(degs, 10))
print(top_10_DE)
view(top_10_DE)

view(degs) #1,510

#Number of upregulated genes
upregulated_genes <- subset(degs, log2FoldChange > 1)
num_upregulated <- nrow(upregulated_genes)
cat("Number of upregulated genes:", num_upregulated, "\n")

#Number of downregulated genes
downregulated_genes <- subset(degs, log2FoldChange < -1)
num_downregulated <- nrow(downregulated_genes)
cat("Number of downregulated genes:", num_downregulated, "\n")




# 2h treatment group

annotation = getBM(attributes=c('ensembl_gene_id', 'chromosome_name',
                                'start_position', 'end_position',
                                'strand', 'gene_biotype', 'external_gene_name',
                                'description'),
                   filters = 'ensembl_gene_id', values = filtered_results2$ensembl_gene_id,
                   mart = ensembl108)
annot_results2 = left_join(filtered_results2, annotation)
annot_results2 = arrange(annot_results2, padj)

degs2 = filter(annot_results2, abs(log2FoldChange) > 1 & padj < 0.05)
top_10_DE2 <- (head(degs2, 10))
print(top_10_DE2)
view(top_10_DE2)

view(degs2) #1,268

#Number of upregulated genes
upregulated_genes2 <- subset(degs2, log2FoldChange > 0)
num_upregulated2 <- nrow(upregulated_genes2)
cat("Number of upregulated genes:", num_upregulated2, "\n")

#Number of downregulated genes
downregulated_genes2 <- subset(degs2, log2FoldChange < 0)
num_downregulated2 <- nrow(downregulated_genes2)
cat("Number of downregulated genes:", num_downregulated2, "\n")





##plot volcano with new colours and x axis limits
# 24h treatment group
volcano_1 <- ggplot(filtered_results, aes(x = log2FoldChange, y = logPVal)) +
  geom_point(aes(colour = combined_significance), alpha = 0.8, size = 2) +
  scale_color_manual(
    values = c(
      'downregulated' = 'cornflowerblue',
      'non-significant' = 'grey',
      'upregulated' = 'indianred2',
      'p value only' = '#E69F00')) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "black") +
  theme_minimal() + xlim(-23, 23) +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)) +
  ggtitle("Differentially expressed genes in Allo24h treatment group") +
  labs(colour = "Differential expression")

volcano_1

#Add gene annotations for top 10 most differentially expressed genes
volcano_24genes <- volcano_1 +
  geom_label_repel(data = top_10_DE,
                   aes(x = log2FoldChange, y = logPVal, label = external_gene_name),
                   box.padding = 0.7, point.padding = 0.5,
                   segment.color = "black", segment.size = 0.5,
                   max.overlaps = Inf) +
  theme(legend.position = "right")
volcano_24genes



#2h treatment group

volcano_2 <- ggplot(filtered_results2, aes(x = log2FoldChange, y = logPVal)) +
  geom_point(aes(colour = combined_significance), alpha = 0.8, size = 2) +
  scale_color_manual(
    values = c(
      'downregulated' = 'cornflowerblue',
      'non-significant' = 'grey',
      'upregulated' = 'indianred2',
      'p value only' = '#E69F00')) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "black") +
  theme_minimal() + xlim(-24, 24) +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)) +
  ggtitle("Differentially expressed genes in Allo2h treatment group") +
  labs(colour = "Differential expression")

volcano_2

#Add gene annotations for top 10 most differentially expressed genes
volcano_2genes <- volcano_2 +
  geom_label_repel(data = top_10_DE2,
                   aes(x = log2FoldChange, y = logPVal, label = external_gene_name),
                   box.padding = 0.7, point.padding = 0.5,
                   segment.color = "black", segment.size = 0.5,
                   max.overlaps = Inf,
                   fill = 'white') +
  theme(legend.position = "right")
volcano_2genes

#####################
