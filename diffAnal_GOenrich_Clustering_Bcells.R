## Load packages ----
library(tidyverse)
library(ggplot2)
library(plotly)
library(magrittr)
library(limma)
library(Glimma)
library(edgeR)
library(FactoMineR)
library(factoextra)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(cluster)
library(mclust)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ReactomePA)
library(enrichplot)
library(igraph)
library(Factoshiny)
library(VennDiagram)
library(pheatmap)



## Functions ----------------------------------------------
discretise <- function(x, t) {
  # Low level function to return a 3 level discretisation.
  if (x <= -t) return(-1)
  else if (x >= t) return(1)
  else return(0)
}

discretiseList <- function(v, thres, fc = 1, ...) {
  # Low level function which discretises a list, either by a 'hard' threshold, or by MAD.
  med <- median(v)
  md <- mad(v)
  #med = mean(v)
  #md = sd(v)
  if (missing(thres)) {
    thres2 = med + fc*md
  } else {
    thres2 = thres
  }
  return(sapply(v, FUN = discretise, t = thres2, ...))
}

discretiseMatrix <- function(df, ...){
  # Function to discretise a matrix.
  ll <- apply(df, 2, FUN = discretiseList, ...)
  return(as.data.frame(ll))
}

plotMatrixBoxplot <- function(df, ...) {
  # Function that plots a violin - jitter plot of a numeric matrix.
  dff <- df %>% rownames_to_column(var = "GeneID") %>% gather(Experiment, LogFC, -GeneID, factor_key = TRUE)
  p <- ggplot(dff, aes(x = Experiment, y = LogFC, color = Experiment)) + geom_violin(trim = FALSE) + geom_jitter(aes(alpha = 0.5), position = position_jitter(0.25)) + stat_summary(fun.data=mean_sdl, geom = "pointrange", color = "dimgrey", size = 1) + stat_summary(fun.y = median, geom = "point", shape = 95, size = 10, color="black") + scale_color_brewer(palette = "Dark2") + theme_linedraw()
  return(p)
}

geneBoxplotCond <- function(matrix, name, experiments, treatments, jit_width = 0.1, point_size = 2, ...){
  # Function to plot expression of individual gene.
  # Experiment: are the different fractions.
  # Treatment: is High or Low glucose.
  ge <- data.frame(t(matrix[name,]));
  ge$exp <- experiments;
  ge$treat <- treatments;
  colnames(ge)[1] <- "TPM";
  p <- ggplot(ge, aes(exp, TPM));
  p + geom_jitter(aes(color = treat), width = jit_width, size = point_size) + ggtitle(name);
}

lowExpression_filter <- function(e, groups, thres = 3, samples = 1, coefVar = 0.5, ...){
  # Function to filter for lowly expressed genes.
  # e      : raw counts data.frame (or cpm, or tpm matrix).
  # groups : factor designating the grouping(s) (conditions, treatments etc.) it MUSt be of equal length to the columns of e and its levels must be the different groups.
  # thres  : the threshold of the *mean* between groups.
  # samples: the minimum number of sample groups that we want threshold 'thres' to be higher.
  # coefVar: The coefficient of variation threshold used to remove "noisy" genes between replicates.
  filteredDF <- data.frame()
  rows <- vector()
  for (i in 1:nrow(e)) {
    row <- as.numeric(e[i,])
    means <- aggregate(row~groups, FUN = mean)$row
    sds <- aggregate(row~groups, FUN = sd)$row
    cvs <- sds/means
    lg <- length(levels(groups))
    if ( (sum(cvs <= coefVar) == lg) & (sum(means >= thres) >= samples) ) { # The first part of this condition checks for the coefficient of variability in ALL groups by grouping sums in the number of samples and the second part is checking for the experiment threshold we specifie.
      filteredDF <- rbind(filteredDF, row)
      rows <- append(rows, rownames(e)[i])
    }
  }
  colnames(filteredDF) <- colnames(e)
  rownames(filteredDF) <- rows
  return(filteredDF)
}


plot_semiSupervised_clust <- function(data, k, method, scale = FALSE, ...){
  # Nicely plots a k-means clustering.
  # Scaling.
  if (scale == TRUE) {
    data <- as.data.frame(scale(data))
  }
  # Calculate the clustering.
  clustMeth <- match.fun(method)
  clustRes <- clustMeth(data, k, ...)
  # Append id and cluster
  dfcall <- cbind(data, id = seq(nrow(data)), cluster = clustRes$cluster)
  # Add idsort, the id number ordered by cluster
  dfcall$idsort <- dfcall$id[order(dfcall$cluster)]
  dfcall$idsort <- order(dfcall$idsort)
  # Generate cluster colours.
  clusterCols <- as.character(sort(clustRes$cluster))
  # Plotting
  heatmap(as.matrix(data)[order(clustRes$cluster),], Rowv = NA, Colv = NA, scale = "none", labRow = NA, cexCol = 1.5, col = my_palette, RowSideColors = clusterCols, ylab = "Genes", main = paste(method, " clustering with ", k, " clusters scale:", as.character(scale)))
  invisible(list(res = clustRes, df = dfcall))
}


plot_unSupervised_clust <- function(data, method, scale = FALSE, ...){
  # Nicely plots a k-means clustering.
  # Scaling.
  if (scale == TRUE) {
    data <- as.data.frame(scale(data))
  }
  # Calculate the clustering.
  clustMeth <- match.fun(method)
  clustRes <- clustMeth(data, ...)
  # Append id and cluster
  dfcall <- cbind(data, id = seq(nrow(data)), cluster = clustRes$classification)
  # Add idsort, the id number ordered by cluster
  dfcall$idsort <- dfcall$id[order(dfcall$cluster)]
  dfcall$idsort <- order(dfcall$idsort)
  # Generate cluster colours.
  clusterCols <- as.character(sort(clustRes$classification))
  # Plotting
  heatmap(as.matrix(data)[order(clustRes$classification),], Rowv = NA, Colv = NA, scale = "none", labRow = NA, cexCol = 1.5, col = my_palette, RowSideColors = clusterCols, ylab = "Genes", main = paste(method, " clustering, scale:", as.character(scale)))
  invisible(list(res = clustRes, df = dfcall))
}




## Preprpocess data DE -------------
# Load the counts table.
countsTableRaw <- read.table("countsTOTALS_CodingGenes.tsv", header = TRUE, sep = "\t")
# Compute the CPM
cpmall <- as.data.frame(cpm(countsTableRaw))

# Divide the Counts Table
countsTableRawT <- countsTableRaw[,1:6]
countsTableRawM <- countsTableRaw[,7:12]
countsTableRawL <- countsTableRaw[,13:18]
countsTableRawH <- countsTableRaw[,19:24]

# Form the proper factors for naming.
groupsformatrix_M <- factor(c("MonoH", "MonoH", "MonoH", "MonoL", "MonoL", "MonoL"),
                            levels = c("MonoH", "MonoL"))
groupsformatrix_L <- factor(c("LightH", "LightH", "LightH", "LightL", "LightL", "LightL"),
                            levels = c("LightH", "LightL"))
groupsformatrix_H <- factor(c("HeavyH", "HeavyH", "HeavyH", "HeavyL", "HeavyL", "HeavyL"),
                            levels = c("HeavyH", "HeavyL"))
groupsformatrix_T <- factor(c("TotalH", "TotalH", "TotalH", "TotalL", "TotalL", "TotalL"),
                            levels = c("TotalH", "TotalL"))

groupsall <- factor(c("Mono", "Mono", "Mono", "Mono", "Mono", "Mono",
                      "Light", "Light", "Light", "Light", "Light", "Light",
                      "Heavy", "Heavy", "Heavy", "Heavy", "Heavy", "Heavy",
                      "Total", "Total", "Total", "Total", "Total", "Total"),
                    levels = c("Mono", "Light", "Heavy", "Total"));

treasAll <- as.factor(c("High", "High", "High", "Low", "Low", "Low"))  ## WHY this thing exists two times we do not know.
treasAll <- as.factor(c(rep(c("High", "High", "High", "Low", "Low", "Low"), 4)));

# Compute the CPM for each division.
cpmtmpM <- cpm(countsTableRawM)
cpmtmpL <- cpm(countsTableRawL)
cpmtmpH <- cpm(countsTableRawH)
cpmtmpT <- cpm(countsTableRawT)

# Remove genes with CPM < 1
cpmtmpM <- cpmtmpM[apply(cpmtmpM > 1, 1, all),]
cpmtmpL <- cpmtmpL[apply(cpmtmpL > 1, 1, all),]
cpmtmpH <- cpmtmpH[apply(cpmtmpH > 1, 1, all),]
cpmtmpT <- cpmtmpT[apply(cpmtmpT > 1, 1, all),]

# Filter the CPM,
cpmallM <- lowExpression_filter(cpmtmpM, groupsformatrix_M, thres = 4, samples = 2)
cpmallL <- lowExpression_filter(cpmtmpL, groupsformatrix_L, thres = 4, samples = 2)
cpmallH <- lowExpression_filter(cpmtmpH, groupsformatrix_H, thres = 4, samples = 2)
cpmallT <- lowExpression_filter(cpmtmpT, groupsformatrix_T, thres = 4, samples = 2)

# Collect all the names of filtred genes.
cpm_Filt_names <- unique(c(rownames(cpmallM), rownames(cpmallL), rownames(cpmallH), rownames(cpmallT)))

# Extract the counts of the filtered genes.
countsTable_H <- countsTableRawH[rownames(cpmallH),]
countsTable_L <- countsTableRawL[rownames(cpmallL),]
countsTable_M <- countsTableRawM[rownames(cpmallM),]
countsTable_T <- countsTableRawT[rownames(cpmallT),]
countsTableAll <- countsTableRaw[cpm_Filt_names,]
cpmAll_Filt <- cpmall[cpm_Filt_names, ]


### Quality Control Plots ------------------------------
## PCA on counts.
res.pca.Counts = PCA(t(countsTableAll), graph = FALSE)
fviz_pca_ind(res.pca.Counts,
             fill.ind = groupsall, col.var = "black", repel = TRUE,
             col.ind = groupsall, # colorer by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "orange"),
             addEllipses = TRUE, # Ellipses de concentration
             legend.title = "Groups",
             title = "PCA on the Polysome profile counts table.")

## PCA on CPMs
res.pca.CPMs = PCA(t(cpmAll_Filt), graph = FALSE)
fviz_pca_ind(res.pca.CPMs,
             fill.ind = groupsall, col.var = "black", repel = TRUE,
             col.ind = groupsall, # colorer by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "orange"),
             addEllipses = TRUE, # Ellipses de concentration
             legend.title = "Groups",
             title = "PCA on the Polysome profile CPMs table.")



## Differential Expression analysis -----------------------
# Form the design matrices
designMat_M <- model.matrix(~0+groupsformatrix_M)
colnames(designMat_M) <- levels(groupsformatrix_M)

designMat_L <- model.matrix(~0+groupsformatrix_L)
colnames(designMat_L) <- levels(groupsformatrix_L)

designMat_H <- model.matrix(~0+groupsformatrix_H)
colnames(designMat_H) <- levels(groupsformatrix_H)

designMat_T <- model.matrix(~0+groupsformatrix_T)
colnames(designMat_T) <- levels(groupsformatrix_T)

# Form the contrast matrices
contr.matrix_M <- makeContrasts(
  Monosomes_HvsL = MonoH - MonoL,
  levels = colnames(designMat_M))

contr.matrix_L <- makeContrasts(
  LightPoly_HvsL = LightH - LightL,
  levels = colnames(designMat_L))

contr.matrix_H <- makeContrasts(
  HeavyPoly_HvsL = HeavyH - HeavyL,
  levels = colnames(designMat_H))

contr.matrix_T <- makeContrasts(
  Total_HvsL = TotalH - TotalL,
  levels = colnames(designMat_T))


### DiffExp Monosomes -------------------------------
vmM <- voom(countsTable_M, designMat_M, plot = TRUE)
fitM <- lmFit(vmM, designMat_M)
vfitM <- contrasts.fit(fitM, contrasts = contr.matrix_M)
efitM <- eBayes(vfitM, robust = TRUE) # Playing with the parameters of ebayes makes no difference.
efM <- decideTests(efitM, p.value = 0.05, lfc = 0.5)
summary(efM)
plotSA(efitM, main = "SA plot for monosomes L/H")
tfitM <- treat(vfitM, lfc = log2(1.1))
ttM <- decideTests(tfitM)
summary(ttM)

plotMD(efitM, column = 1, status = efM[,1], main = colnames(efitM)[1],ylim = c( -1.5, 1.5))


### DiffExp Light ribosomes -------------------------------
vmL <- voom(countsTable_L, designMat_L, plot = TRUE)
fitL <- lmFit(vmL, designMat_L)
vfitL <- contrasts.fit(fitL, contrasts = contr.matrix_L)
efitL <- eBayes(vfitL, robust = TRUE) # Playing with the parameters of ebayes makes no difference.
efL <- decideTests(efitL, p.value = 0.05, lfc = 0.5)
summary(efL)
plotSA(efitL, main = "SA plot for light polysomes L/H")
tfitL <- treat(vfitL, lfc = log2(1.1))
ttL <- decideTests(tfitL)
summary(ttL)

plotMD(efitL, column = 1, status = efL[,1], main = colnames(efitL)[1], ylim = c( -1.5, 1.5))


### DiffExp Heavy ribosomes -------------------------------
vmH <- voom(countsTable_H, designMat_H, plot = TRUE)
fitH <- lmFit(vmH, designMat_H)
vfitH <- contrasts.fit(fitH, contrasts = contr.matrix_H)
efitH <- eBayes(vfitH, robust = TRUE) # Playing with the parameters of ebayes makes no difference.
efH <- decideTests(efitH,p.value = 0.05, lfc = 0.5)
summary(efH)
plotSA(efitH, main = "SA plot for heavy polysomes L/H")
tfitH <- treat(vfitH, lfc = log2(1.1))
ttH <- decideTests(tfitH)
summary(ttH)

plotMD(efitH, column = 1, status = efH[,1], main = colnames(efitH)[1],ylim = c( -1.5, 1.5))


### DiffExp Total RNA -------------------------------
vmT <- voom(countsTable_T, designMat_T, plot = TRUE)
fitT <- lmFit(vmT, designMat_T)
vfitT <- contrasts.fit(fitT, contrasts = contr.matrix_T)
efitT <- eBayes(vfitT, robust = TRUE) # Playing with the parameters of ebayes makes no difference.
efT <- decideTests(efitT,p.value = 0.05, lfc = 0.5)
summary(efT)
plotSA(efitT, main = "SA plot for RNA total L/H")
tfitT <- treat(vfitT, lfc = log2(1.15))
ttT <- decideTests(tfitT)
summary(ttT)

plotMD(efitT, column = 1, status = efT[,1], main = colnames(efitT)[1],ylim = c( -1.5, 1.5))

### Toptables ------------
DEH <- topTable(efitH, coef = 1, p.value = 0.05, lfc = 0.5, number = Inf)
DEL <- topTable(efitL, coef = 1, p.value = 0.05, lfc = 0.5, number = Inf)
DEM <- topTable(efitM, coef = 1, p.value = 0.05, lfc = 0.5, number = Inf)
DET <- topTable(efitT, coef = 1, p.value = 0.05, lfc = 0.5, number = Inf)

degs <- union(rownames(DEL), rownames(DEH))
write(degs, "degs_09042019_geneNames.txt")

### Venn diagrams -----------
# Intersection between Light, Heavy, Total
vennALL <- venn.diagram(list(Heavy = rownames(DEH), Light = rownames(DEL), Total = rownames(DET)), NULL, fill = c("darkorange1", "deepskyblue3", "darkolivegreen4"), alpha = c(0.5, 0.5, 0.5), cex = 3)



#### Preprocess Translation data -------------------------
tpmPPglu <- read.table("polysomeProfile_TPM_proteinCoding.csv", header = TRUE, sep = ";")

row_NON_zero <- apply(tpmPPglu, 1, function(row) all(row != 0))
tpmPPgluClean <- tpmPPglu[row_NON_zero,]

# Then we collect the genes and remove the total.
tpmDEGs <- tpmPPglu[degs,]
tpmDEGs <- tpmDEGs[,1:18]

# Calculate means of each group
monoHdeg <- rowMeans(tpmDEGs[,1:3])
monoLdeg <- rowMeans(tpmDEGs[,4:6])
lightHdeg <- rowMeans(tpmDEGs[,7:9])
lightLdeg <- rowMeans(tpmDEGs[,10:12])
heavyHdeg <- rowMeans(tpmDEGs[,13:15])
heavyLdeg <- rowMeans(tpmDEGs[,16:18])

monoHall <- rowMeans(tpmPPgluClean[,1:3])
monoLall <- rowMeans(tpmPPgluClean[,4:6])
lightHall <- rowMeans(tpmPPgluClean[,7:9])
lightLall <- rowMeans(tpmPPgluClean[,10:12])
heavyHall <- rowMeans(tpmPPgluClean[,13:15])
heavyLall <- rowMeans(tpmPPgluClean[,16:18])

# Calculate the ratios.
ratioHall <- log2(heavyHall/heavyLall)
ratioLall <- log2(lightHall/lightLall)
ratioMall <- log2(monoHall/monoLall)

ratioHdeg <- log2(heavyHdeg/heavyLdeg)
ratioLdeg <- log2(lightHdeg/lightLdeg)
ratioMdeg <- log2(monoHdeg/monoLdeg)

# Put the logRatios all in a dataframe.
logRatiosDEG <- data.frame("Mono" = ratioMdeg, "Light" = ratioLdeg, "Heavy" = ratioHdeg,  row.names = rownames(tpmDEGs))
logRatiosALL <- data.frame("Mono" = ratioMall, "Light" = ratioLall, "Heavy" = ratioHall,  row.names = rownames(tpmPPgluClean))



## Enrichments ------------------------------------------
# Generate the gene list.
geneListENS <- logRatiosDEG[["Heavy"]]
names(geneListENS) <- degs
geneListENS <- sort(geneListENS, decreasing = TRUE)
genesENS <- degs
# Generate a gene list with gene symbols too
geneConv <- bitr(degs, fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db)
rownames(geneConv) <- geneConv[["ENSEMBL"]]
geneConv$ENSEMBL <- NULL
genesSYMB <- vector()
for (name in names(geneListENS)) {
  nSymb <- geneConv[name,]
  genesSYMB <- append(genesSYMB, nSymb)
}
geneListSYMB <- geneListENS
names(geneListSYMB) <- genesSYMB

# Create the "universe" gene set (this is a set of all the "detected" genes)
row_NON_zero <- apply(tpmPPglu, 1, function(row) all(row != 0))
tpmPPgluClean <- tpmPPglu[row_NON_zero,]
univPPglu <- rownames(tpmPPgluClean)
geneConv2 <- bitr(univPPglu, fromType = "ENSEMBL", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db)
univPPglu2 <- geneConv2[["SYMBOL"]]

## MSIGDB enrichment
m_df <- msigdbr(species = "Homo sapiens")
m_df$gs_id <- m_df$gene_symbol # BIG TRICK TO SWAP THE COLUMN NAMES!!!!!!

# Gene enrichment
esigDEGs <- enricher(genesSYMB, universe = univPPglu2, TERM2GENE = m_df, minGSSize = 50, maxGSSize = 1000)  # One can play arounf with the set sizes to obtain something meaningfull.
barplot(esigDEGs, showCategory = 50, title = "BarPlot DEGs MsiGDB enrichment")
dotplot(esigDEGs, showCategory = 50, title = "DotPlot DEGs MsiGDB enrichment")

# Gene set enrichment
esigsDEGs <- GSEA(geneListSYMB, minGSSize = 20, TERM2GENE = m_df)
dotplot(esigsDEGs, showCategory = 50, title = "DotPlot MsiGDB DEGS GSEA")


### GO enrichments ------------
# GO group enrichment
ggoDEGs_MF2 <- groupGO(gene = genesENS, OrgDb = org.Hs.eg.db, ont = "MF", level = 2, keyType = "ENSEMBL", readable = TRUE)
barplot(ggoDEGs_MF2, showCategory = 30,  title = "GroupGO DEGs MF_2")
ggoDEGs_MF4 <- groupGO(gene = genesENS, OrgDb = org.Hs.eg.db, ont = "MF", level = 4, keyType = "ENSEMBL", readable = TRUE)
barplot(ggoDEGs_MF4, showCategory = 60,  title = "GroupGO DEGs MF_4")

ggoDEGs_BP2 <- groupGO(gene = genesENS, OrgDb = org.Hs.eg.db, ont = "BP", level = 2, keyType = "ENSEMBL", readable = TRUE)
barplot(ggoDEGs_BP2, showCategory = 30,  title = "GroupGO DEGs BP_2")
ggoDEGs_BP4 <- groupGO(gene = genesENS, OrgDb = org.Hs.eg.db, ont = "BP", level = 4, keyType = "ENSEMBL", readable = TRUE)
barplot(ggoDEGs_BP4, showCategory = 100,  title = "GroupGO DEGs BP_4")

# GO enrichment
egoDEGs_MF <- enrichGO(gene = genesENS, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
dotplot(egoDEGs_MF, title = "GO enrichment DEGs MF")
egoDEGs_BP <- enrichGO(gene = genesENS, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
dotplot(egoDEGs_MF, title = "GO enrichment DEGs BP")
# and one for ALL
egoDEGs_ALL <- enrichGO(gene = genesENS, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", pool = TRUE, readable = TRUE)
dotplot(egoDEGs_ALL, title = "ALL GO enrichment DEGs")

# GO GSEA enrichment
egogsDEGs_MF <- gseGO(geneList = geneListENS, OrgDb = org.Hs.eg.db, ont = "MF", nPerm = 1000, pvalueCutoff = 0.05, keyType = "ENSEMBL")
dotplot(egogsDEGs_MF, title = "GSEA GO MF DEGs")
egogsDEGs_BP <- gseGO(geneList = geneListENS, OrgDb = org.Hs.eg.db, ont = "BP", nPerm = 1000, pvalueCutoff = 0.05, keyType = "ENSEMBL")
dotplot(egogsDEGs_BP, title = "GSEA GO BP DEGs")
egogsDEGs_ALL <- gseGO(geneList = geneListENS, OrgDb = org.Hs.eg.db, ont = "ALL", nPerm = 1000, pvalueCutoff = 0.05, keyType = "ENSEMBL")
dotplot(egogsDEGs_ALL, title = "GSEA GO ALL DEGs")

#TODO include the WikiPathways enrichments also!

## Perform KEEG pathway enrichment.
# Convert ENSEMBL IDs to ENTREZ
genesENTREZ <- as.character(mapIds(org.Hs.eg.db, genesENS, "ENTREZID", "ENSEMBL"))
# Enrich KEGG pathways
ekegDEGs <- enrichKEGG(gene = genesENTREZ, organism = "hsa", pvalueCutoff = 0.05)
barplot(ekegDEGs, title = "DEGs KEGG enrichment")

# Enrich KEGG modules
ekegMDGEs <- enrichMKEGG(gene = genesENTREZ, organism = "hsa", pvalueCutoff = 0.05)
barplot(ekegMDGEs, title = "DEGs KEGG modules enrichment")

# Enrich REACTOME Pathways
ekePDEGs <- enrichPathway(gene = genesENTREZ, organism = "human", pvalueCutoff = 0.05)
barplot(ekePDEGs, showCategory = 30, title = "DEGs REACTOME Pathways enrichment")


### Enrichment Visualisation ------------------------------

## Category Network (CNET) plots (perhaps the most usefull!)
cnetplot(egoDEGs_MF, foldChange = geneListENS, colorEdge = TRUE) + ggtitle("CNETplot GOenrich DEGs MF")
cnetplot(egoDEGs_BP, foldChange = geneListENS, colorEdge = TRUE) + ggtitle("CNETplot GOenrich DEGs BP")
cnetplot(egoDEGs_ALL, foldChange = geneListENS, colorEdge = TRUE) + ggtitle("CNETplot GOenrich DEGs ALL")

cnetplot(egogsDEGs_MF, foldChange = geneListENS, colorEdge = TRUE) + ggtitle("CNETplot GOgsea DEGs MF")
cnetplot(egogsDEGs_BP, foldChange = geneListENS, colorEdge = TRUE) + ggtitle("CNETplot GOgsea DEGs BP")
cnetplot(egogsDEGs_ALL, foldChange = geneListENS, colorEdge = TRUE) + ggtitle("CNETplot GOgsea DEGs ALL")



#### Clustering -------------------------------------------
my_palette <- brewer.pal(n = 11, name = "RdYlGn")
c
### Hierarchical Clustering -------------------------------
# Clustering.
hc_S <- hclust(dist(logRatiosDEG), method = "single")
# define clusters (hard thresold)
hc_S_Cls <- cutree(hc_S, h = max(hc_S$height/4))
# Colour vector for clusters side bar.
myCols_hc_S <- rainbow(length(unique(hc_S_Cls)))
myClusts_hc_S <- myCols_hc_S[hc_S_Cls]
heatmap.2(as.matrix(logRatiosDEG), main = "DEGs logRatio H/L HClust Single",  Rowv = as.dendrogram(hc_S), Colv = FALSE, dendrogram = "row", col = my_palette, cexCol = 1.5, cexRow = 0.5, key.title = NA, keysize = 0.8, key.xlab = NA, ylab = "Genes", RowSideColors = myClusts_hc_S)

# Independent clustering just for comparison.
heatmap.2(as.matrix(logRatiosDEG), main = "DEGs logRatio H/L HClust Average",  dendrogram = "row", Colv = FALSE, col = my_palette, cexCol = 1.5, cexRow = 0.5, key.title = NA, keysize = 0.8, key.xlab = NA, ylab = "Genes", hclustfun = function(x) hclust(x, method = "average"))

heatmap.2(as.matrix(logRatiosDEG), main = "DEGs logRatio H/L HClust Ward",  dendrogram = "row", Colv = FALSE, col = my_palette, cexCol = 1.5, cexRow = 0.5, key.title = NA, keysize = 0.8, key.xlab = NA, ylab = "Genes", hclustfun = function(x) hclust(x, method = "ward.D2"))

heatmap.2(as.matrix(logRatiosDEG), main = "DEGs logRatio H/L HClust Single",  dendrogram = "row", Colv = FALSE, col = my_palette, cexCol = 1.5, cexRow = 0.5, key.title = NA, keysize = 0.8, key.xlab = NA, ylab = "Genes", hclustfun = function(x) hclust(x, method = "single"))

heatmap.2(as.matrix(logRatiosDEG), main = "DEGs logRatio H/L HClust Complete",  dendrogram = "row", Colv = FALSE, col = my_palette, cexCol = 1.5, cexRow = 0.5, key.title = NA, keysize = 0.8, key.xlab = NA, ylab = "Genes", hclustfun = function(x) hclust(x, method = "complete"))


## Visualise results from different clustering algorithms.
# Hierarchical clustering
clustH.res <- hclust(dist(logRatiosDEG), method = "ward.D2")
fviz_cluster(clustH.res, data = logRatiosDEG, ellipse.type = "convex") + theme_minimal()

# K-means
kmeans4.res <- kmeans(scale(logRatiosDEG), 4)
kmeans5.res <- kmeans(scale(logRatiosDEG), 5)
plot_semiSupervised_clust(logRatiosDEG, 4, "kmeans")
fviz_cluster(kmeans5.res, data = logRatiosDEG, ellipse.type = "convex") + theme_minimal()

# PAM
pam4.res <- pam(logRatiosDEG, 4)
pam5.res <- pam(logRatiosDEG, 5)
plot_semiSupervised_clust(logRatiosDEG, 4, "pam")
plot_semiSupervised_clust(logRatiosDEG, 5, "pam")
fviz_cluster(pam5.res, data = logRatiosDEG, ellipse.type = "convex") + theme_minimal()

# MClust
mclust.res <- Mclust(logRatiosDEG)
fviz_cluster(mclust.res, data = logRatiosDEG, ellipse.type = "convex") + theme_minimal()

### We choose the Mclust method!
clustRes <- (plot_unSupervised_clust(logRatiosDEG, "Mclust"))

## 3D interactive plot of the MClust result.




## RNA Features analysis --------
# Prepare the features data frame.
clusterGeneIDs <- clustRes$df["cluster"]
featuresDF <- read.table("rnaFeat/degs_09042019_ENSEMBL.tab", header = TRUE, sep = ";")
# Create a slice with only the numeric values of the data frame.
featDF <- featuresDF[,c(1:13)]

# Add the clustering column
featDF["Cluster"] <- 0
for (r in row.names(featDF)) {
  geneID <- as.character(featDF[r,]$ensembl_gene_id);
  featDF[r,]$Cluster <- clusterGeneIDs[geneID,];
}

# Add the translation column.
featDF["Translation"] <- "a"

# Characterise the clusters based on translation behaviour.
#! ATTENTION !# This is ONLY for the specific clustering result!
for (r in row.names(featDF)) {
  if (featDF[r,]$Cluster %in% c(1,6)) featDF[r,]$Translation <- "Up"
  if (featDF[r,]$Cluster %in% c(4,5)) featDF[r,]$Translation <- "Down"
  if (featDF[r,]$Cluster %in% c(2,3)) featDF[r,]$Translation <- "Inter"
}

# Plot the Coding length boxplots
ggplot(featDF, aes(x = Cluster, y = coding_len, fill = Translation, group = Cluster)) +
  coord_cartesian(ylim = c(0, 4000)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.95, notch = TRUE) +
  theme(legend.position = "topleft") +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("Coding Length") +
  xlab("Cluster") +
  ggtitle("Coding length distribution among clusters") +
  theme_bw()

# Plot the 5UTR length boxplots
ggplot(featDF, aes(x = Cluster, y = len_5pUTR, fill = Translation, group = Cluster)) +
  coord_cartesian(ylim = c(0, 1000)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.95, notch = TRUE) +
  theme(legend.position = "topleft") +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("5'UTR Length") +
  xlab("Cluster") +
  ggtitle("5'UTR length distribution among clusters") +
  theme_bw()

# Plot the GC 5'UTR boxplots
ggplot(featDF, aes(x = Cluster, y = GC_5pUTR, fill = Translation, group = Cluster)) +
  coord_cartesian(ylim = c(35, 100)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.95, notch = TRUE) +
  theme(legend.position = "topleft") +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("5'UTR GC") +
  xlab("Cluster") +
  ggtitle("GC 5'UTR distribution among clusters") +
  theme_bw()

# Plot the 5UTR MFE boxplots
ggplot(featDF, aes(x = Cluster, y = MFE_5pUTR, fill = Translation, group = Cluster)) +
  coord_cartesian(ylim = c(-400, 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.95, notch = TRUE) +
  theme(legend.position = "topleft") +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("5'UTR MFE") +
  xlab("Cluster") +
  ggtitle("MFE 5'UTR distribution among clusters") +
  theme_bw()

# Plot the 5UTR MFE_BP boxplots
ggplot(featDF, aes(x = Cluster, y = MfeBP_5pUTR, fill = Translation, group = Cluster)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.95, notch = TRUE) +
  theme(legend.position = "topleft") +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("5'UTR Mfe_Bp") +
  xlab("Cluster") +
  ggtitle("MFE pre Bp 5'UTR distribution among clusters") +
  theme_bw()

# Plot the 3UTR lengths boxplots
ggplot(featDF, aes(x = Cluster, y = len_3pUTR, fill = Translation, group = Cluster)) +
  coord_cartesian(ylim = c(0, 4000)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.95, notch = TRUE) +
  theme(legend.position = "topleft") +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("3'UTR length") +
  xlab("Cluster") +
  ggtitle("3'UTR length distribution among clusters") +
  theme_bw()

# Plot the 3UTR GC lengths boxplots #TODO FIX THE COLUMN NAME
ggplot(featDF, aes(x = Cluster, y = len_3pUTR.1, fill = Translation, group = Cluster)) +
  #coord_cartesian(ylim = c(0, 4000)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.95, notch = TRUE) +
  theme(legend.position = "topleft") +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("3'UTR GC") +
  xlab("Cluster") +
  ggtitle("3'UTR GC distribution among clusters") +
  theme_bw()

# Plot the 3UTR MFE lengths boxplots
ggplot(featDF, aes(x = Cluster, y = MFE_3pUTR, fill = Translation, group = Cluster)) +
  coord_cartesian(ylim = c(-1500, 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.95, notch = TRUE) +
  theme(legend.position = "topleft") +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("3'UTR MFE") +
  xlab("Cluster") +
  ggtitle("3'UTR MFE distribution among clusters") +
  theme_bw()

# Plot the 3UTR MFE_BP lengths boxplots
ggplot(featDF, aes(x = Cluster, y = MfeBP_3pUTR, fill = Translation, group = Cluster)) +
  #coord_cartesian(ylim = c(-1500, 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.95, notch = TRUE) +
  theme(legend.position = "topleft") +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("3'UTR MFE_BP") +
  xlab("Cluster") +
  ggtitle("3'UTR MFE per BP distribution among clusters") +
  theme_bw()

# Plot the TOP local score boxplots
ggplot(featDF, aes(x = Cluster, y = TOP_localScore, fill = Translation, group = Cluster)) +
  #coord_cartesian(ylim = c(-1500, 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.95, notch = TRUE) +
  theme(legend.position = "topleft") +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("TOP local score") +
  xlab("Cluster") +
  ggtitle("TOP local score distribution among clusters") +
  theme_bw()

# Plot the CAI index boxplots
ggplot(featDF, aes(x = Cluster, y = CAI, fill = Translation, group = Cluster)) +
  #coord_cartesian(ylim = c(-1500, 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.95, notch = TRUE) +
  theme(legend.position = "topleft") +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("CAI index") +
  xlab("Cluster") +
  ggtitle("CAI index distribution among clusters") +
  theme_bw()
