## Load packages ------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggstatsplot)
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
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(msigdbr)
library(ReactomePA)
library(enrichplot)
library(igraph)
library(Factoshiny)
library(pheatmap)
library(VennDiagram)
library(eulerr)
library(heatmap3)
library(ggpubr)
library(GGally)

## Functions ----------------------------------------------
# Source functions from the .r file.
source("/home/costas/devel/u-paris/diffAnalysis_GOenrich_Clustering_Bcells/diffAnal_GOenrich_Clustering_Bcells_FunctionSource.r")


## Preprocess data DE ------------------------------------
# Load the counts table.
countsTableRaw <- read.table("data/countsTOTALS_CodingGenes.tsv", header = TRUE, sep = "\t")

# Compute the CPM
cpmall <- as.data.frame(cpm(countsTableRaw))

# Divide the Counts Table
countsTableRawT <- countsTableRaw[,1:6]
countsTableRawM <- countsTableRaw[,7:12]
countsTableRawL <- countsTableRaw[,13:18]
countsTableRawH <- countsTableRaw[,19:24]

# Form the proper factors for naming.
groupsformatrix_M <- factor(c("MonoH", "MonoH", "MonoH", "MonoL", "MonoL", "MonoL"), levels = c("MonoH", "MonoL"))
groupsformatrix_L <- factor(c("LightH", "LightH", "LightH", "LightL", "LightL", "LightL"), levels = c("LightH", "LightL"))
groupsformatrix_H <- factor(c("HeavyH", "HeavyH", "HeavyH", "HeavyL", "HeavyL", "HeavyL"), levels = c("HeavyH", "HeavyL"))
groupsformatrix_T <- factor(c("TotalH", "TotalH", "TotalH", "TotalL", "TotalL", "TotalL"), levels = c("TotalH", "TotalL"))

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
cpmallM <- filter_low_expression(cpmtmpM, groupsformatrix_M, thres = 4, samples = 2)
cpmallL <- filter_low_expression(cpmtmpL, groupsformatrix_L, thres = 4, samples = 2)
cpmallH <- filter_low_expression(cpmtmpH, groupsformatrix_H, thres = 4, samples = 2)
cpmallT <- filter_low_expression(cpmtmpT, groupsformatrix_T, thres = 4, samples = 2)

# Collect all the names of filtered genes.
cpm_Filt_names <- unique(c(rownames(cpmallM), rownames(cpmallL), rownames(cpmallH), rownames(cpmallT)))

# Extract the counts of the filtered genes.
countsTable_H <- countsTableRawH[rownames(cpmallH),]
countsTable_L <- countsTableRawL[rownames(cpmallL),]
countsTable_M <- countsTableRawM[rownames(cpmallM),]
countsTable_T <- countsTableRawT[rownames(cpmallT),]
countsTableAll <- countsTableRaw[cpm_Filt_names,]
cpmAll_Filt <- cpmall[cpm_Filt_names, ]


### Quality Control Plots ---------------------------------
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

# Cluster transcriptome CPMs.
hcT <- hclust(dist(t(cpmallT)), method = "ward.D2")
clus2 = cutree(hcT, 2)
colCl <- c("green", "red")
plot(as.dendrogram(hcT))
fviz_dend(hcT, cex = 1.5, main = "", k =2, k_colors = c("#00AFBB", "#FC4E07"), color_labels_by_k = TRUE, xlab = "Total mRNA samples", ylab = "", sub = "", ggtheme = theme_ggstatsplot())


## Differential Expression analyses -----------------------
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


### DiffExp Monosomes -------------------------------------
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

plotMD(efitM, column = 1, status = efM[,1], main = colnames(efitM)[1], ylim = c( -1.5, 1.5))
# For the manuscript figure.
par(mar = c(7,6,3,0))
plotMD(efitM, column = 1, status = efM[,1], ylim = c( -1.5, 1.5),  xlab = "Average log occupancy", ylab= "log Fold-Change", main = "Monosome occupancy on High vs. Low glucose", hl.col = c("lightgreen", "red3"), cex.lab = 2, cex.axis = 2, cex = 2, cex.main = 2, legend = FALSE)


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
# For the manuscript figure.
par(mar = c(7,6,3,0))
plotMD(efitL, column = 1, status = efL[,1], ylim = c( -1.5, 1.5), xlab = "Average log occupancy", ylab= "log Fold-Change", main = "Light-polysomes occupancy on High vs. Low glucose", hl.cex = 2,  hl.col = c("lightgreen", "red3"), cex.lab = 2, cex.axis = 2, cex = 2, cex.main = 1.75, legend = FALSE)

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
# For the manuscript figure.
par(mar = c(7,6,3,0))
plotMD(efitH, column = 1, status = efH[,1], ylim = c( -1.5, 1.5), xlab = "Average log occupancy", ylab= "log Fold-Change", main = "Heavy-polysomes occupancy on High vs. Low glucose", hl.cex = 2,  hl.col = c("lightgreen", "red3"), cex.lab = 2, cex.axis = 2, cex = 2, cex.main = 1.75, legend = FALSE)


### DiffExp Total RNA -------------------------------------
vmT <- voom(countsTable_T, designMat_T, plot = TRUE)
fitT <- lmFit(vmT, designMat_T)
vfitT <- contrasts.fit(fitT, contrasts = contr.matrix_T)
efitT <- eBayes(vfitT, robust = TRUE) # Playing with the parameters of ebayes makes no difference.
efT <- decideTests(efitT, p.value = 0.05, lfc = 0.5)
summary(efT)
efitTtab <- topTable(efitT, p.value = 0.05, lfc = 0.5, number = "ALL")
efitTtab$GeneName <- as.vector(mapIds(org.Hs.eg.db, keys = rownames(efitTtab), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
write.table(efitTtab, file = "degs_RNAtotal_201911.tab", quote = FALSE, sep = "\t")
plotSA(efitT, main = "SA plot for RNA total L/H")
tfitT <- treat(vfitT, lfc = log2(1.1))
ttT <- decideTests(tfitT, p.value = 0.05)
summary(ttT)

plotMD(efitT, column = 1, status = efT[,1], main = colnames(efitT)[1], ylim = c( -1.1, 1.1))
plotMD(tfitT, column = 1, status = ttT[,1], main = colnames(tfitT)[1], ylim = c( -1.1, 1.1))
# For the manuscript figure.
par(mar = c(7,6,3,0))
plotMD(efitT, column = 1, status = efT[,1], ylim = c( -1, 1), xlab = "Average log expression", ylab= "log Fold-Change", main = "Total RNA abundance on High vs. Low glucose", hl.cex = 2,  hl.col = c("lightgreen", "red3"), cex.lab = 2, cex.axis = 2, cex = 2, cex.main = 1.75, legend = FALSE)


### Toptables ---------------------------------------------
DEH <- topTable(efitH, coef = 1, p.value = 0.05, lfc = 0.5, number = Inf)
DEL <- topTable(efitL, coef = 1, p.value = 0.05, lfc = 0.5, number = Inf)
DEM <- topTable(efitM, coef = 1, p.value = 0.05, lfc = 0.5, number = Inf)
DET <- topTable(efitT, coef = 1, p.value = 0.05, lfc = 0.5, number = Inf)
DEHup <- subset(DEH, DEH$logFC > 0)
DEHdown <- subset(DEH, DEH$logFC < 0)
DELup <- subset(DEL, DEL$logFC > 0)
DELdown <- subset(DEL, DEL$logFC < 0)
DETup <- subset(DET, DET$logFC > 0)
DETdown <- subset(DET, DET$logFC < 0)


degs <- union(rownames(DEL), rownames(DEH))
write(degs, "degs_17092020_geneNames.txt")


### Venn diagrams -----------------------------------------
# Intersection between Light, Heavy, Total
# DEPRECATED vennALL <- venn.diagram(list(Heavy = rownames(DEH), Light = rownames(DEL), Total = rownames(DET)), NULL, fill = c("darkorange1", "deepskyblue3", "darkolivegreen4"), alpha = c(0.5, 0.5, 0.5), cex = 3) DEPRECATED #

# Venn diagrams with eulerr package
polyVenn <- list("Heavy_UP" = rownames(DEHup), "Light_UP" = rownames(DELup), "Total_UP" = rownames(DETup), "Heavy_DOWN" = rownames(DEHdown), "Light_DOWN" = rownames(DELdown), "Total_DOWN" = rownames(DETdown))

plot(euler(polyVenn, shape = "ellipse"), quantities = TRUE)

# VEnn diagram figure 1 paper.


#### Process Translation data --------------------------
tpmPPglu <- read.table("data/polysomeProfile_TPM_proteinCoding.csv", header = TRUE, sep = ";")

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


## Enrichments --------------------------------------------
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


### GO enrichments ----------------------------------------
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
dotplot(egoDEGs_BP, title = "GO enrichment DEGs BP")
# and one for ALL
egoDEGs_ALL <- enrichGO(gene = genesENS, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", pool = TRUE, readable = TRUE)
dotplot(egoDEGs_ALL, title = "ALL GO enrichment DEGs")

# GO GSEA enrichment
egogsDEGs_MF <- gseGO(geneList = geneListSYMB, OrgDb = org.Hs.eg.db, ont = "MF", nPerm = 1000, pvalueCutoff = 0.05, keyType = "SYMBOL")
dotplot(egogsDEGs_MF, title = "GSEA GO MF DEGs")
egogsDEGs_BP <- gseGO(geneList = geneListSYMB, OrgDb = org.Hs.eg.db, ont = "BP", nPerm = 1000, pvalueCutoff = 0.05, keyType = "SYMBOL")
dotplot(egogsDEGs_BP, title = "GSEA GO BP DEGs")
egogsDEGs_ALL <- gseGO(geneList = geneListSYMB, OrgDb = org.Hs.eg.db, ont = "ALL", nPerm = 1000, pvalueCutoff = 0.05, keyType = "SYMBOL")
dotplot(egogsDEGs_ALL, title = "GSEA GO ALL DEGs")

#TODO include the WikiPathways enrichments also!

## Perform KEEG pathway enrichment.
# Convert ENSEMBL IDs to ENTREZ
genesENTREZ <- as.character(mapIds(org.Hs.eg.db, genesENS, "ENTREZID", "ENSEMBL"))
# Enrich KEGG pathways
ekegDEGs <- enrichKEGG(gene = genesENTREZ, organism = "hsa", pvalueCutoff = 0.05)
barplot(ekegDEGs, title = "DEGs KEGG enrichment")  # Only one "ribosome"

# Enrich KEGG modules
ekegMDGEs <- enrichMKEGG(gene = genesENTREZ, organism = "hsa", pvalueCutoff = 0.05)
barplot(ekegMDGEs, title = "DEGs KEGG modules enrichment")  # Only one "ribosome"

# Enrich REACTOME Pathways
ekePDEGs <- enrichPathway(gene = genesENTREZ, organism = "human", pvalueCutoff = 0.05)
barplot(ekePDEGs, showCategory = 30, title = "DEGs REACTOME Pathways enrichment")
dotplot(ekePDEGs, showCategory = 20, font.size = 18)  #, title = "DEGs REACTOME Pathways enrichment")


### Enrichment Visualisation ------------------------------
# Category Network (CNET) plots (perhaps the most usefull!)
cnetplot(egoDEGs_MF, foldChange = geneListENS, colorEdge = TRUE, showCategory = 10) + ggtitle("CNETplot GOenrich DEGs MF")
cnetplot(egoDEGs_BP, foldChange = geneListENS, colorEdge = TRUE, showCategory = 10) + ggtitle("CNETplot GOenrich DEGs BP")
cnetplot(egoDEGs_ALL, foldChange = geneListENS, colorEdge = TRUE, showCategory = 10) + ggtitle("CNETplot GOenrich DEGs ALL")
# Paper figure 2A
egoDEGs_MF_paper <- egoDEGs_MF
egoDEGs_MF_paper@result <- egoDEGs_MF@result[-c(2, 9),]  # Remove two redundant categories
cnetplot(egoDEGs_MF_paper, foldChange = geneListENS, colorEdge = TRUE, showCategory = 8,cex_category =2, cex_label_category = 1.5, cex_gene =0.75, shadowtext= 'category')
# Paper figure 2B
egoDEGs_BP_paper <- egoDEGs_BP
egoDEGs_BP_paper@result <- egoDEGs_BP@result[-c(1,2,5,6),]  # remove 4 redundant categories
cnetplot(egoDEGs_BP_paper, foldChange = geneListENS, colorEdge = TRUE, showCategory = 6, cex_category = 1.5, cex_label_category = 1.4, cex_gene = 0.5, cex_label_gene = 1, shadowtext= 'category')

cnetplot(egogsDEGs_MF, foldChange = geneListSYMB, colorEdge = TRUE, showCategory = 10) + ggtitle("CNETplot GOgsea DEGs MF")
cnetplot(egogsDEGs_BP, foldChange = geneListSYMB, colorEdge = TRUE, showCategory = 10) + ggtitle("CNETplot GOgsea DEGs BP")
cnetplot(egogsDEGs_ALL, foldChange = geneListSYMB, colorEdge = TRUE, showCategory = 10) + ggtitle("CNETplot GOgsea DEGs ALL")

cnetplot(ekePDEGs, foldChange = geneListSYMB, colorEdge = TRUE, symbol = "ENSEMBL") + ggtitle("CNETplot GOgsea DEGs ALL")


#### Clustering -------------------------------------------
my_palette <- brewer.pal(n = 11, name = "RdYlGn")

### Hierarchical Clustering -------------------------------
# Clustering.
hc_S <- hclust(dist(logRatiosDEG), method = "single")
# define clusters (hard threshold)
hc_S_Cls <- cutree(hc_S, h = max(hc_S$height/4))
# Colour vector for clusters side bar.
myCols_hc_S <- rainbow(length(unique(hc_S_Cls)))mclust
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

# MClust!
mclust.res <- Mclust(logRatiosDEG)
fviz_cluster(mclust.res, data = logRatiosDEG, ellipse.type = "convex") + theme_minimal()

### We choose the MClust method!
plot_unSupervised_clust(logRatiosDEG, "Mclust")

# The paper figure.
par(mar = c(0,0,0,0))
heatmap.2(as.matrix(logRatiosDEG)[order(mclust.res$classification),], Colv = NA, Rowv = NA, labRow = NA, cexCol = 4, col = my_palette, RowSideColors = brewer.pal(n = 6, name = "Dark2")[as.factor(as.character(sort(mclust.res$classification)))], trace = "none", lwid = c(0.8, 3), lhei = c(0.7, 3), margins = c(9,0), srtCol = 45, key.ylab = "Density", key.xlab = "logOdds", key.title = "")
#mtext("Genes", side = 4, cex = 2.5, padj = -1.5)
#mtext(colnames(logRatiosDEG), cex = 2, side = 1, srt = 45, padj = -5)
legend("bottomleft", legend=c(1:6), fill = brewer.pal(6, "Dark2"), border = brewer.pal(6, "Dark2"), cex = 2.5, bty = "n", y.intersp = 1.7, title = "Cluster")

# TODO 3D interactive plot of the MClust result.


## RNA Features analysis ----------------------------------
# Prepare the features data frame.q
clusterGeneIDs <- clustRes$df["cluster"]
featuresDF <- read.table("rnaFeat/201905/degs_02052019_ENSEMBL.tab", header = TRUE, sep = "\t")
# Create a slice with only the numeric values of the data frame.
featDF <- featuresDF[,c(1:15)]

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

# Features correlations


## BOXPLOTS ---------

# Coding length boxplots
clCD <- ggplot(featDF[featDF$coding_len <= 4000,], aes(x = Cluster, y = coding_len, fill = Cluster)) + theme_pubr() + scale_fill_brewer(palette="Dark2") +
  #ggtitle("Coding length distributions of the 6 MClust clusters")
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.1, height = 0)) +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  ylab("Coding Length") + xlab("Clusters")
  #stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
  #stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(6, 2), c(6, 3), c(6, 4), c(6, 5), c(2, 4), c(3, 5)), method = "t.test", tip.length = 0, hide.ns = TRUE, size = 2.5)

clCDa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(), ylab = "Coding Length",
  data = featDF[featDF$coding_len <= 4000,],  # Select the ones less than 4000.
  x = Cluster, y = coding_len,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "Coding length distributions of the 6 MClust clusters")

trCDa <- ggbetweenstats( # Group translation
  ggtheme = theme_ggstatsplot(), ylab = "Coding Length", xlab = "",
  data = featDF[featDF$coding_len <= 4000,],
  x = Translation, y = coding_len,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "Coding length distributions of the 3 translation behaviours")

# GC content boxplots
clGC <- ggplot(featDF, aes(x = Cluster, y = GC, fill = Cluster)) + theme_pubr() +
  scale_fill_brewer(palette="Dark2") +
  coord_cartesian(ylim = c(30, 75)) +
  #ggtitle("Transcript GC content distributions of the 6 Mclust clusters") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.1, height = 0)) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("GC") + xlab("Clusters")
  #stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
  #stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(6, 2), c(6, 3), c(6, 4), c(6, 5), c(2, 4), c(3, 5)), method = "t.test", tip.length = 0, hide.ns = TRUE, size = 2.5)

clGCa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featDF, ylab = "GC",
  x = Cluster, y = GC,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 1, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "GC content distributions of the 6 MClust clusters")

trGCa <- ggbetweenstats( # Group translation
  ggtheme = theme_ggstatsplot(),
  data = featDF[featDF$coding_len <= 4000,], ylab = "GC",
  x = Translation, y = GC,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 1, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "GC content distributions of the 3 translation behaviours")

# 5'UTR length boxplots
cl5pLEN <- ggplot(featDF[featDF$len_5pUTR <= 800 & featDF$len_5pUTR != 0,], aes(x = Cluster, y = len_5pUTR, fill = Cluster)) + theme_pubr() +
  scale_fill_brewer(palette="Dark2") +
  #ggtitle("5'UTR length distributions of the 6 Mclust cluster") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.1, height = 0)) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("5'UTR length") + xlab("Clusters")
#stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
#stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(6, 2), c(6, 3), c(6, 4), c(6, 5), c(2, 4), c(3, 5)), method = "t.test", tip.length = 0, hide.ns = TRUE, size = 2.5)

cl5pLENa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featDF[featDF$len_5pUTR <= 1000,], ylab = "5'UTR length",
  x = Cluster, y = len_5pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "5'UTR length distributions of the 6 MClust clusters")

tr5pLENa <- ggbetweenstats( # Group translation
  ggtheme = theme_ggstatsplot(),
  data = featDF[featDF$len_5pUTR <= 1000,], ylab = "5'UTR length",
  x = Translation, y = len_5pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "5'UTR length distributions of the 3 translation behaviours")

# 5'UTR GC boxplots
cl5pGC <- ggplot(featDF[featDF$len_5pUTR != 0,], aes(x = Cluster, y = GC_5pUTR, fill = Cluster)) + theme_pubr() +
  scale_fill_brewer(palette="Dark2") +
  #ggtitle("5'UTR GC distributions of the 6 Mclust cluster") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.1, height = 0)) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("5'UTR GC") + xlab("Clusters")
#stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
#stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(6, 2), c(6, 3), c(6, 4), c(6, 5), c(2, 4), c(3, 5)), method = "t.test", tip.length = 0, hide.ns = TRUE, size = 2.5)

cl5pGCa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featDF[featDF$len_5pUTR > 1,], ylab = "GC 5pUTR",
  x = Cluster, y = GC_5pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 1, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "5'UTR GC distributions of the 6 MClust clusters")

tr5pGCa <-  ggbetweenstats( # Group translation
  ggtheme = theme_ggstatsplot(),
  data = featDF[featDF$GC_5pUTR >= 25,], ylab = "GC 5pUTR",
  x = Translation, y = GC_5pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 1, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "5'UTR GC distributions of the 3 translation behaviours")

# 5'UTR MFE boxplots
cl5pMFE <- ggplot(featDF[featDF$MFE_5pUTR > -400 & featDF$len_5pUTR != 0,], aes(x = Cluster, y = MFE_5pUTR, fill = Cluster)) + theme_pubr() +
  scale_fill_brewer(palette="Dark2") +
  #coord_cartesian(ylim = c(-400, 0)) +
  #ggtitle("5'UTR MFE distributions of the 6 Mclust cluster") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.1, height = 0)) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("5'UTR MFE") #+ xlab("Clusters")
#stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
#stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(6, 2), c(6, 3), c(6, 4), c(6, 5), c(2, 4), c(3, 5)), method = "t.test", tip.length = 0, hide.ns = TRUE, size = 2.5)

cl5pMFEa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featDF[featDF$MFE_5pUTR >= -400,], ylab = "5'UTR MFE",
  x = Cluster, y = MFE_5pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "5'UTR MFE distributions of the 6 MClust clusters")

tr5pMFEa <- ggbetweenstats( # Group translation
  ggtheme = theme_ggstatsplot(),
  data = featDF[featDF$MFE_5pUTR >= -400,], ylab = "5'UTR MFE",
  x = Translation, y = MFE_5pUTR,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "MFE 5'UTR distributions of the 3 translation behaviours")

# 5UTR MFE_BP boxplots
cl5pMfeBP <- ggplot(featDF[featDF$len_5pUTR > 1,], aes(x = Cluster, y = MfeBP_5pUTR, fill = Cluster)) + theme_pubr() +
  scale_fill_brewer(palette="Dark2") +
  #ggtitle("5'UTR MFE per BP distributions of the 6 Mclust cluster") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.1, height = 0)) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("5'UTR MFE") #+ xlab("Clusters")
#stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
#stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(6, 2), c(6, 3), c(6, 4), c(6, 5), c(2, 4), c(3, 5)), method = "t.test", tip.length = 0, hide.ns = TRUE, size = 2.5)

cl5pMfeBPa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featDF[featDF$MfeBP_5pUTR != 0,], ylab = "5'UTR MFE/BP",
  x = Cluster, y = MfeBP_5pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 2, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "MFE per BP 5'UTR distributions of the 6 MClust clusters")

tr5pMfeBPa <- ggbetweenstats( # Group translation
  ggtheme = theme_ggstatsplot(),
  data = featDF[featDF$MfeBP_5pUTR != 0,], ylab = "5'UTR MFE/BP",
  x = Translation, y = MfeBP_5pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 2, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "MFE per BP 5'UTR distributions of the 3 translation behaviours")

# 3UTR lengths boxplots
cl3pLEN <- ggplot(featDF[featDF$len_3pUTR <= 4000 & featDF$len_3pUTR !=0,], aes(x = Cluster, y = len_3pUTR, fill = Cluster)) + theme_pubr() +
  scale_fill_brewer(palette="Dark2") +
  #ggtitle("3'UTR length distributions of the 6 Mclust cluster") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.1, height = 0)) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("3'UTR length") #+ xlab("Clusters")
#stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
#stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(6, 2), c(6, 3), c(6, 4), c(6, 5), c(2, 4), c(3, 5)), method = "t.test", tip.length = 0, hide.ns = TRUE, size = 2.5)

cl3pLENa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(), ylab = "3'UTR length",
  data = featDF[featDF$len_3pUTR != 0 & featDF$len_3pUTR <= 4000,],
  x = Cluster, y = len_3pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "3'UTR length distributions of the 6 MClust clusters")

tr3pLENa <- ggbetweenstats( # Group translation
  ggtheme = theme_ggstatsplot(), ylab = "3'UTR length",
  data = featDF[featDF$len_3pUTR != 0 & featDF$len_3pUTR <= 4000,],
  x = Translation, y = len_3pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "3'UTR length distributions of the 3 translation behaviours")

# 3UTR GC boxplots
cl3pGC <- ggplot(featDF[featDF$len_3pUTR !=0,], aes(x = Cluster, y = GC_3pUTR, fill = Cluster)) + theme_pubr() +
  scale_fill_brewer(palette="Dark2") +
  #ggtitle("3'UTR GC distributions of the 6 Mclust cluster") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.1, height = 0)) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("3'UTR GC") #+ xlab("Clusters")
#stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
#stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(6, 2), c(6, 3), c(6, 4), c(6, 5), c(2, 4), c(3, 5)), method = "t.test", tip.length = 0, hide.ns = TRUE, size = 2.5)

cl3pGCa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(), ylab = "3'UTR GC",
  data = featDF[featDF$len_3pUTR != 0,],
  x = Cluster, y = GC_3pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 1, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "GC 3'UTR distributions of the 6 MClust clusters")

tr3pGCa <- ggbetweenstats( # Group translation
  ggtheme = theme_ggstatsplot(), ylab = "3'UTR GC",
  data = featDF[featDF$len_3pUTR != 0,],
  x = Translation, y = GC_3pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 1, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "GC 3'UTR distribution of the 3 translation behaviours")

# Plot the 3UTR MFE boxplots
cl3pMFE <- ggplot(featDF[featDF$MFE_3pUTR >= -1500 & featDF$len_3pUTR !=0,], aes(x = Cluster, y = MFE_3pUTR, fill = Cluster)) + theme_pubr() +
  scale_fill_brewer(palette="Dark2") +
  coord_cartesian(ylim = c(-1500, 0)) +
  #ggtitle("3'UTR MFE distributions of the 6 Mclust cluster") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.1, height = 0)) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("3'UTR MFE") + xlab("Clusters")
#stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
#stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(6, 2), c(6, 3), c(6, 4), c(6, 5), c(2, 4), c(3, 5)), method = "t.test", tip.length = 0, hide.ns = TRUE, size = 2.5)

cl3pMFEa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(), ylab = "3'UTR MFE",
  data = featDF[featDF$MFE_3pUTR >= -1500 & featDF$MFE_3pUTR != 0,],
  x = Cluster, y = MFE_3pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "MFE 3'UTR distributions of the 6 MClust clusters")

tr3pMFEa <- ggbetweenstats( # Group translation
  ggtheme = theme_ggstatsplot(), ylab = "3'UTR MFE",
  data = featDF[featDF$MFE_3pUTR >= -1500 & featDF$MFE_3pUTR != 0,],
  x = Translation, y = MFE_3pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "MFE 3'UTR distributions of the 3 translation behaviours")

# Plot the 3UTR MFE_BP lengths boxplots
cl3pMfeBP <- ggplot(featDF[featDF$MFE_3pUTR >= -1500 & featDF$len_3pUTR !=0,], aes(x = Cluster, y = MfeBP_3pUTR, fill = Cluster)) + theme_pubr() +
  scale_fill_brewer(palette="Dark2") +
  #ggtitle("3'UTR MFE per BP distributions of the 6 Mclust cluster") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.1, height = 0)) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("3'UTR MFE per BP") + xlab("Clusters")
#stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
#stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(6, 2), c(6, 3), c(6, 4), c(6, 5), c(2, 4), c(3, 5)), method = "t.test", tip.length = 0, hide.ns = TRUE, size = 2.5)

cl3pMfeBPa <- ggbetweenstats( # Group clusters
  ggtheme =theme_pubr(), ylab = "3'UTR MFE/BP",
  data = featDF[featDF$MFE_3pUTR >= -1500 & featDF$MFE_3pUTR != 0,],
  x = Cluster, y = MfeBP_3pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 2, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "MFE per BP 3'UTR distributions of the 6 MClust clusters")

tr3pMfeBPa <- ggbetweenstats( # Group translation
  ggtheme = theme_ggstatsplot(), ylab = "3'UTR MFE/BP",
  data = featDF[featDF$MFE_3pUTR >= -1500 & featDF$MFE_3pUTR != 0,],
  x = Translation, y = MfeBP_3pUTR,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 2, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "MFE per BP 3'UTR distributions of the 3 translation behaviours")

# Plot the TOP local score boxplots
clTOP <- ggplot(featDF[featDF$len_5pUTR !=0,], aes(x = Cluster, y = TOP_localScore, fill = Cluster)) + theme_pubr() +
  scale_fill_brewer(palette="Dark2") +
  #ggtitle("TOP local score distributions of the 6 Mclust cluster") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.1, height = 0.05)) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("TOP local score") + xlab("Cluster")
#stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
#stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(6, 2), c(6, 3), c(6, 4), c(6, 5), c(2, 4), c(3, 5)), method = "t.test", tip.length = 0, hide.ns = TRUE, size = 2.5)

clTOPa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featDF, ylab = "TOP local score",
  x = Cluster, y = TOP_localScore,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 1, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "TOP local score distributions of the 6 MClust clusters")

trTOPa <- ggbetweenstats( # Group translation
  ggtheme = theme_ggstatsplot(),
  data = featDF, ylab = "TOP local score",
  x = Translation, y = TOP_localScore,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 1, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "TOP local score distributions of the 3 translation behaviours")

# Plot the CAI index boxplots
clCAI <- ggplot(featDF, aes(x = Cluster, y = CAI, fill = Cluster)) + theme_pubr() +
  scale_fill_brewer(palette="Dark2") +
  coord_cartesian(ylim = c(0.65, 0.9)) +
  #ggtitle("CAI distributions of the 6 Mclust cluster") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.1, height = 0.05)) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  scale_x_discrete(limits = c("1","2","3","4","5","6")) +
  ylab("CAI") + xlab("Cluster")
#stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
#stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(6, 2), c(6, 3), c(6, 4), c(6, 5), c(2, 4), c(3, 5)), method = "t.test", tip.length = 0, hide.ns = TRUE, size = 2.5)

clCAIa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featDF, ylab = "CAI",
  x = Cluster, y = CAI,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 2, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "CAI distributions of the 6 MClust clusters")

trCAIa <- ggbetweenstats( # Group translation
  ggtheme = theme_ggstatsplot(),
  data = featDF, ylab = "CAI",
  x = Translation, y = CAI,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 2, mean.point.args = list(size = 1.5, color = "darkred"),
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE,  # Set FALSE for the manuscript figure.
  title = "CAI distribution of the 3 translation behaviours")

# Produce combined publication figures.
figureBoxClustP <- ggarrange(clCD, clGC, cl5pLEN, cl5pGC, cl5pMFE, cl5pMfeBP,
                             cl3pLEN, cl3pGC, cl3pMFE, cl3pMfeBP, clTOP,
                             clCAI,
                             labels = c("A", "B", "C", "D", "E" ,"F", "G", "H",
                                        "I", "J", "K", "L"),
                             ncol = 2, nrow = 6)

analysisBoxClustI <- ggarrange(clCDa, clGCa, cl5pLENa, cl5pGCa, cl5pMFEa,
                              cl5pMfeBPa,
                              labels = c("A", "B", "C", "D", "E" ,"F"),
                              ncol = 2, nrow = 3)

analysisBoxClustII <- ggarrange(cl3pLENa, cl3pGCa, cl3pMFEa,
                                cl3pMfeBPa, clTOPa, clCAIa,
                                labels = c("G", "H", "I", "J", "K", "L"),
                                ncol =2, nrow = 3)

analysisBoxTranslI <- ggarrange(trCDa, trGCa, tr5pLENa, tr5pGCa, tr5pMFEa,
                              tr5pMfeBPa,
                              labels = c("A", "B", "C", "D", "E" ,"F"),
                              ncol = 2, nrow = 3)

analysisBoxTranslII <- ggarrange(tr3pLENa, tr3pGCa, tr3pMFEa,
                                tr3pMfeBPa, trTOPa, trCAIa,
                                labels = c("G", "H", "I", "J", "K", "L"),
                                ncol = 2, nrow = 3)


## Cluster Enrichements ------------------------------------
# Retrieve the gene names for each cluster from the features data frame.
genesClust1ENS <- as.vector(subset(featDF, featDF$Cluster == 1)$ensembl_gene_id)
genesClust2ENS <- as.vector(subset(featDF, featDF$Cluster == 2)$ensembl_gene_id)
genesClust3ENS <- as.vector(subset(featDF, featDF$Cluster == 3)$ensembl_gene_id)
genesClust4ENS <- as.vector(subset(featDF, featDF$Cluster == 4)$ensembl_gene_id)
genesClust5ENS <- as.vector(subset(featDF, featDF$Cluster == 5)$ensembl_gene_id)
genesClust6ENS <- as.vector(subset(featDF, featDF$Cluster == 6)$ensembl_gene_id)
genesTranslUP <- as.vector(subset(featDF, featDF$Translation == "Up")$ensembl_gene_id)
genesTranslDOWN <- as.vector(subset(featDF, featDF$Translation == "Down")$ensembl_gene_id)
genesTranslINTER <- as.vector(subset(featDF, featDF$Translation == "Inter")$ensembl_gene_id)

# GO enrichments
#Cluster1
egoClust1_MF <- enrichGO(gene = genesClust1ENS, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoClust1_MF, title = "Cluster1 GO enrichment MF", showCategory = 20)
egoClust1_BP <- enrichGO(gene = genesClust1ENS, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01,   minGSSize = 5, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoClust1_BP, showCategory = 30, font.size = 16) #title = "Cluster1 GO BP enrichment",
egoClust1_ALL <- enrichGO(gene = genesClust1ENS, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", pool = TRUE, readable = TRUE)
barplot(egoClust1_ALL, title = "Cluster1 ALL GO enrichment", showCategory = 20)

#Cluster6
egoClust6_MF <- enrichGO(gene = genesClust6ENS, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoClust6_MF, title = "Cluster6 GO enrichment MF", showCategory = 20)
egoClust6_BP <- enrichGO(gene = genesClust6ENS, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, maxGSSize = 1000, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoClust6_BP, showCategory = 30, font.size = 18)  # title = "Cluster6 GO BP enrichment",
egoClust6_ALL <- enrichGO(gene = genesClust6ENS, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", pool = TRUE, readable = TRUE)
barplot(egoClust6_ALL, title = "Cluster6 ALL GO enrichment", showCategory = 20)

#Cluster4
egoClust4_MF <- enrichGO(gene = genesClust4ENS, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoClust4_MF, title = "Cluster4 GO enrichment MF", showCategory = 20)
egoClust4_BP <- enrichGO(gene = genesClust4ENS, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, minGSSize = 6, maxGSSize = 1000, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoClust4_BP, showCategory = 30, font.size = 18)  # title = "Cluster4 GO BP enrichment",
egoClust4_ALL <- enrichGO(gene = genesClust4ENS, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", pool = TRUE, readable = TRUE)
barplot(egoClust4_ALL, title = "Cluster4 ALL GO enrichment", showCategory = 20)

#Cluster5
egoClust5_MF <- enrichGO(gene = genesClust5ENS, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoClust5_MF, title = "Cluster5 GO enrichment MF", showCategory = 20)
egoClust5_BP <- enrichGO(gene = genesClust5ENS, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, minGSSize = 5, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoClust5_BP, showCategory = 20, font.size = 14)  # title = "Cluster5 GO BP enrichment",
egoClust5_ALL <- enrichGO(gene = genesClust5ENS, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", pool = TRUE, readable = TRUE)
barplot(egoClust5_ALL, title = "Cluster5 ALL GO enrichment", showCategory = 20)

#Cluster2
egoClust2_MF <- enrichGO(gene = genesClust2ENS, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoClust2_MF, title = "Cluster2 GO enrichment MF", showCategory = 20)
egoClust2_BP <- enrichGO(gene = genesClust2ENS, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.1, minGSSize = 5, maxGSSize = 1000, keyType = "ENSEMBL", readable = TRUE)
barplot(egoClust2_BP, title = "Cluster2 GO BP enrichment", showCategory = 20)
egoClust2_ALL <- enrichGO(gene = genesClust2ENS, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", pool = TRUE, readable = TRUE)
barplot(egoClust2_ALL, title = "Cluster2 ALL GO enrichment", showCategory = 20)

#Cluster3
egoClust3_MF <- enrichGO(gene = genesClust3ENS, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoClust3_MF, title = "Cluster3 GO enrichment MF", showCategory = 20)
egoClust3_BP <- enrichGO(gene = genesClust3ENS, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.1, minGSSize = 5, maxGSSize = 1000, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoClust3_BP, title = "Cluster3 GO BP enrichment", showCategory = 20)
egoClust3_ALL <- enrichGO(gene = genesClust3ENS, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", pool = TRUE, readable = TRUE)
barplot(egoClust3_ALL, title = "Cluster3 ALL GO enrichment", showCategory = 20)

#Translation UP
egoTranslUP_MF <- enrichGO(gene = genesTranslUP, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoTranslUP_MF, title = "TranslUP GO enrichment MF", showCategory = 20)
egoTranslUP_BP <- enrichGO(gene = genesTranslUP, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoTranslUP_BP, title = "TranslUP GO BP enrichment", showCategory = 20)
egoTranslUP_ALL <- enrichGO(gene = genesTranslUP, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", pool = TRUE, readable = TRUE)
barplot(egoTranslUP_ALL, title = "TranslUP ALL GO enrichment", showCategory = 20)

#Translation Inter
egoTranslInter_MF <- enrichGO(gene = genesTranslINTER, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoTranslInter_MF, title = "TranslInter GO enrichment MF", showCategory = 20)
egoTranslInter_BP <- enrichGO(gene = genesTranslINTER, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoTranslInter_BP, title = "TranslInter GO BP enrichment", showCategory = 20)
egoTranslInter_ALL <- enrichGO(gene = genesTranslINTER, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", pool = TRUE, readable = TRUE)
barplot(egoTranslInter_ALL, title = "TranslInter ALL GO enrichment", showCategory = 20)

#Translation DOWN
egoTranslDOWN_MF <- enrichGO(gene = genesTranslDOWN, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoTranslDOWN_MF, title = "TranslDOWN GO enrichment MF", showCategory = 20)
egoTranslDOWN_BP <- enrichGO(gene = genesTranslDOWN, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", readable = TRUE)
barplot(egoTranslDOWN_BP, title = "TranslDOWN GO BP enrichment", showCategory = 20)
egoTranslDOWN_ALL <- enrichGO(gene = genesTranslDOWN, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = univPPglu, keyType = "ENSEMBL", pool = TRUE, readable = TRUE)
barplot(egoTranslDOWN_ALL, title = "TranslDOWN ALL GO enrichment", showCategory = 20)


## Translation Ratio analyses -----------------------------
# Prepare the datasets.
dd <- read.table("data/polysomeProfile_TPM_proteinCoding.csv", header = TRUE, sep = ";")
rs <- apply(dd, 1, function(row) all(row != 0))
dd <- dd[rs,]
# Filter also for the lowly expressed genes from the diff. analysis before.
dd <- dd[cpm_Filt_names,]
# Get the means of each condition.
monoHavg <- rowMeans(dd[,1:3])
monoLavg <- rowMeans(dd[,4:6])
lightHavg <- rowMeans(dd[,7:9])
lightLavg <- rowMeans(dd[,10:12])
heavyHavg <- rowMeans(dd[,13:15])
heavyLavg <- rowMeans(dd[,16:18])
totalHavg <- rowMeans(dd[,19:21])
totalLavg <- rowMeans(dd[,22:24])
tpmAvgAll <- data.frame("MonoL" = monoLavg, "MonoH" = monoHavg, "LightL" = lightLavg, "LightH" = lightHavg, "HeavyL" = heavyLavg, "HeavyH" = heavyHavg, "TotalL" = totalLavg, "TotalH" = totalHavg)
tpmTrRatAll <- data.frame("MonoL" = monoLavg/totalLavg, "MonoH" = monoHavg/totalHavg, "LightL" = lightLavg/totalLavg, "LightH" = lightHavg/totalHavg, "HeavyL" = heavyLavg/totalLavg, "HeavyH" = heavyHavg/totalHavg)
tpmTrRatAll.m <- reshape2::melt(as.matrix(tpmTrRatAll), id.vars = NULL)
colnames(tpmTrRatAll.m) <- c("GeneID", "Experiment", "TranslRatio")

# First exploratory plot of translation efficiency.
ggplot(tpmTrRatAll.m, aes(x = Experiment, y = TranslRatio, color = Experiment)) + geom_violin(trim = TRUE) + geom_boxplot(width = 0.1, outlier.alpha = 0.2) + theme_minimal() + labs(x = "Experiment", y = "Transl. Ratio") + ylim(0, 7)

## Read all the flatten and pre-processed file.
#tpmTrEffAll.f <- read.table("trans_ratio/tpm_transEff_All_flat.csv", header = TRUE, sep = ",")
# Plot it
#ggplot(tpmTrEffAll.f, aes(x = Experiment, y = TransEfficiency)) + geom_violin(aes(fill = Fraction)) + geom_boxplot(width = 0.1, outlier.alpha = 0.2) + theme_minimal() + labs(x = "Experiment", y = "Trans. Efficiency") + ylim(0, 7)

# Translation ratio, low high glucose.
tpmTrRat <- data.frame("HighGlu" = tpmTrRatAll$LightH + tpmTrRatAll$HeavyH, "LowGlu" = tpmTrRatAll$LightL + tpmTrRatAll$HeavyL)
rownames(tpmTrRat) <- rownames(tpmTrRatAll)
tpmTrRat$GeneID <- rownames(tpmTrRat)

# Translation efficiency difference.
diffTranslRat <- data.frame(difTransRat = (tpmTrRat$HighGlu - tpmTrRat$LowGlu), row.names = rownames(tpmTrRat))

# Keep the 200 most differentially translated gene names in a file.
diffTransl_genes <- rownames(head(diffTranslRat[with(diffTranslRat, order(-difTransRat)),, drop = FALSE], n = 200))
write(diffTransl_genes, file = "most_diffTrans_HiLowGlu.txt")

# Generate a data frame of the average of trans-eff in low and high glucose.
tpmTrRat <- data.frame("HighGlu" = (tpmTrRatAll$LightH + tpmTrRatAll$HeavyH)/2, "LowGlu" = (tpmTrRatAll$LightL + tpmTrRatAll$HeavyL)/2)
rownames(tpmTrRat) <- rownames(tpmTrRatAll)
tpmTrRat.m <- reshape2::melt(as.matrix(tpmTrRat), id.vars = NULL)
colnames(tpmTrRat.m) <- c("GeneID", "Treatment", "TranslRat")

# Have a look of the different translation efficiencies in low and high.
ggplot(tpmTrRat.m, aes(x = Treatment, y = TranslRat)) +
  geom_violin(aes(fill = Treatment), draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5, trim = TRUE) +
  geom_boxplot(width = 0.2, outlier.alpha = 0.2, notch = TRUE) +
  labs(x = "Experiment", y = "Transl. Ratio") + ylim(0, 10)

# Keep the most translated gene names in high glucose in a file
low_TranslRat <- rownames(head(tpmTrRat[with(tpmTrRat, order(-LowGlu)),,], n = 200))
write(low_TranslRat, file = "most_Transl_LowGlu.txt")
high_TranslRat <- rownames(head(tpmTrRat[with(tpmTrRat, order(-HighGlu)),,], n = 200))
write(high_TranslRat, file = "most_Transl_HighGlu.txt")

# Compare all the most differential translated, the most efficiently translated in high and most efficiently translated in low genes.
translRatios <- list(TranslDiff = diffTransl_genes, TranslRatio_Hi = high_TranslRat, TranslRatio_Low = low_TranslRat)
plot(euler(translRatios, shape = "ellipse"), quantities = TRUE)

# Generate the venn object.
vG <- construct(translRatios)

# Read the features for the translation ratio genes.
featuresTransl <- read.table("trans_ratio/translRatio_Features_201906/transl_Ratios_geneNames.tab", header = TRUE, sep = ";")

# Create a slice with only the numeric values of the data frame.
featTransl <- featuresTransl[,c(1:14)]

# Add the clustering column
featTransl["Cluster"] <- NaN

for (i in 1:nrow(featTransl)) {
  if (featTransl[i, "ensembl_gene_id"] %in% overlap(vG)) {
    featTransl[i, "Cluster"] = "All"
  }
  else if (featTransl[i, "ensembl_gene_id"] %in% overlap(vG, c("TranslDiff", "TranslRatio_Hi"))) {
    featTransl[i, "Cluster"] = "HiANDDiff"
  }
  else if (featTransl[i, "ensembl_gene_id"] %in% overlap(vG, c("TranslRatio_Low", "TranslRatio_Hi"))) {
    featTransl[i, "Cluster"] = "LowANDHigh"
  }
  else if (featTransl[i, "ensembl_gene_id"] %in% high_TranslRat) {
    featTransl[i, "Cluster"] = "Hi_Glu"
  }
  else if (featTransl[i, "ensembl_gene_id"] %in% low_TranslRat) {
    featTransl[i, "Cluster"] = "Low_Glu"
  }
  else if (featTransl[i, "ensembl_gene_id"] %in% diffTransl_genes) {
    featTransl[i, "Cluster"] = "Transl_Diff"
  }
}

# Make the plots.
# Plot the Coding length boxplots
ggplot(featTransl, aes(x = Cluster, y = coding_len, fill = Cluster, group = Cluster)) +
  #coord_cartesian(ylim = c(0, 7000)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  geom_jitter(aes(col = Cluster), position = position_jitter(width = .2, height = 0)) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
  theme(legend.position = "topleft") +
  scale_x_discrete(labels = c("All", "Hi Gluc", "Diff + High", "Low Gluc", "Low + High", "Transl Diff")) +
  ylab("Coding Length") +
  xlab("Transl. Groups") +
  ggtitle("Coding length distribution of the 6 Transl. Ratio groups") +
  theme_bw()

# Plot GC content boxplots
ggplot(featTransl, aes(x = Cluster, y = GC, fill = Cluster, group = Cluster)) +
  coord_cartesian(ylim = c(30, 80)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  geom_jitter(aes(col = Cluster), position = position_jitter(width = .2, height = 0)) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "jitter") +
  theme(legend.position = "topleft") +
  scale_x_discrete(labels = c("All", "Hi Gluc", "Diff + High", "Low Gluc", "Low + High", "Transl Diff")) +
  ylab("GC") +
  xlab("Transl. Groups") +
  ggtitle("GC content distribution of the 6 Transl. Ratio groups") +
  theme_bw()

# Plot 5'UTR length boxplots
ggplot(featTransl, aes(x = Cluster, y = X5pUTR_len, fill = Cluster, group = Cluster)) +
  #coord_cartesian(ylim = c(30, 80)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  geom_jitter(aes(col = Cluster), position = position_jitter(width = .2, height = 0)) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "jitter") +
  theme(legend.position = "topleft") +
  scale_x_discrete(labels = c("All", "Hi Gluc", "Diff + High", "Low Gluc", "Low + High", "Transl Diff")) +
  ylab("5'UTR lengths") +
  xlab("Transl. Groups") +
  ggtitle("5'UTR length distribution of the 6 Transl. Ratio groups") +
  theme_bw()

# Plot 5'UTR MFE boxplots
ggplot(featTransl, aes(x = Cluster, y = X5pUTR_MFE, fill = Cluster, group = Cluster)) +
  #coord_cartesian(ylim = c(30, 80)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  geom_jitter(aes(col = Cluster), position = position_jitter(width = .2, height = 0)) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "jitter") +
  theme(legend.position = "topleft") +
  scale_x_discrete(labels = c("All", "Hi Gluc", "Diff + High", "Low Gluc", "Low + High", "Transl Diff")) +
  ylab("5'UTR MFE") +
  xlab("Transl. Groups") +
  ggtitle("5'UTR MFE distribution of the 6 Transl. Ratio groups") +
  theme_bw()

# Plot 5'UTR MFE boxplots
ggplot(featTransl, aes(x = Cluster, y = X5pUTR_MfeBP, fill = Cluster, group = Cluster)) +
  coord_cartesian(ylim = c(0.01, -0.7)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  geom_jitter(aes(col = Cluster), position = position_jitter(width = .2, height = 0)) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "stack") +
  theme(legend.position = "topleft") +
  scale_x_discrete(labels = c("All", "Hi Gluc", "Diff + High", "Low Gluc", "Low + High", "Transl Diff")) +
  ylab("5'UTR MFE/BP") +
  xlab("Transl. Groups") +
  ggtitle("5'UTR MFE/BP distribution of the 6 Transl. Ratio groups") +
  theme_bw()

# Plot 3'UTR length boxplots
ggplot(featTransl, aes(x = Cluster, y = X3pUTR_len, fill = Cluster, group = Cluster)) +
  coord_cartesian(ylim = c(0, 10000)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  geom_jitter(aes(col = Cluster), position = position_jitter(width = .2, height = 0)) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill") +
  theme(legend.position = "topleft") +
  scale_x_discrete(labels = c("All", "Hi Gluc", "Diff + High", "Low Gluc", "Low + High", "Transl Diff")) +
  ylab("3'UTR length") +
  xlab("Transl. Groups") +
  ggtitle("3'UTR length distribution of the 6 Transl. Ratio groups") +
  theme_bw()

# Plot 3'UTR GC boxplots
ggplot(featTransl, aes(x = Cluster, y = X3pUTR_GC, fill = Cluster, group = Cluster)) +
  #coord_cartesian(ylim = c(30, 80)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  geom_jitter(aes(col = Cluster), position = position_jitter(width = .2, height = 0)) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "jitter") +
  theme(legend.position = "topleft") +
  scale_x_discrete(labels = c("All", "Hi Gluc", "Diff + High", "Low Gluc", "Low + High", "Transl Diff")) +
  ylab("3'UTR GC") +
  xlab("Transl. Groups") +
  ggtitle("3'UTR GC distribution of the 6 Transl. Ratio Groups") +
  theme_bw()

# Plot 3'UTR MFE boxplots
ggplot(featTransl, aes(x = Cluster, y = X3pUTR_MFE, fill = Cluster, group = Cluster)) +
  coord_cartesian(ylim = c(0, -2500)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  geom_jitter(aes(col = Cluster), position = position_jitter(width = .2, height = 0)) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "jitter") +
  theme(legend.position = "topleft") +
  scale_x_discrete(labels = c("All", "Hi Gluc", "Diff + High", "Low Gluc", "Low + High", "Transl Diff")) +
  ylab("3'UTR MFE") +
  xlab("Transl. Groups") +
  ggtitle("3'UTR MFE distribution of the 6 Transl. Ratio Groups") +
  theme_bw()

# Plot 3'UTR MFE/BP boxplots
ggplot(featTransl, aes(x = Cluster, y = X3pUTR_MfeBP, fill = Cluster, group = Cluster)) +
  coord_cartesian(ylim = c(0, -0.6)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  geom_jitter(aes(col = Cluster), position = position_jitter(width = .2, height = 0)) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "stack") +
  theme(legend.position = "topleft") +
  scale_x_discrete(labels = c("All", "Hi Gluc", "Diff + High", "Low Gluc", "Low + High", "Transl Diff")) +
  ylab("3'UTR MFE/BP") +
  xlab("Transl. Groups") +
  ggtitle("3'UTR MFE/BP distribution of the 6 Transl. Ratio Groups") +
  theme_bw()

# Plot TOP local score boxplots
ggplot(featTransl, aes(x = Cluster, y = TOP_localScore, fill = Cluster, group = Cluster)) +
  #coord_cartesian(ylim = c(0, -0.6)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  geom_jitter(aes(col = Cluster), position = position_jitter(width = .2, height = 0)) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "dodge") +
    theme(legend.position = "topleft") +
  scale_x_discrete(labels = c("All", "Hi Gluc", "Diff + High", "Low Gluc", "Low + High", "Transl Diff")) +
  ylab("TOP local score") +
  xlab("Transl. Groups") +
  ggtitle("TOP local score distribution of the 6 Transl. Ratio Groups") +
  theme_bw()

# Plot CAI boxplots
ggplot(featTransl, aes(x = Cluster, y = CAI, fill = Cluster, group = Cluster)) +
  #coord_cartesian(ylim = c(0, -0.6)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  geom_jitter(aes(col = Cluster), position = position_jitter(width = .2, height = 0)) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "dodge") +
  theme(legend.position = "topleft") +
  scale_x_discrete(labels = c("All", "Hi Gluc", "Diff + High", "Low Gluc", "Low + High", "Transl Diff")) +
  ylab("CAI") +
  xlab("Transl. Groups") +
  ggtitle("CAI distribution of the 6 Transl. Ratio Groups") +
  theme_bw()


## Analysis of Translation differences --------------------
# Apply a more stringent filter to avoid the lowly translated genes.
# re-read the data
ddT <- read.table("data/polysomeProfile_TPM_proteinCoding.csv", header = TRUE, sep = ";")
# We do not calculate the monosomes fraction.
ddT <- ddT[,7:24]
# Generate the grouping factors
groupsFactor_Trasnl <- factor(c("LightH", "LightH", "LightH", "LightL", "LightL", "LightL",
                                "HeavyH", "HeavyH", "HeavyH", "HeavyL", "HeavyL", "HeavyL",
                                "TotalH", "TotalH", "TotalH", "TotalL", "TotalL", "TotalL"),
                                levels = c("LightH", "LightL", "HeavyH", "HeavyL", "TotalH", "TotalL"));

# Strict filter. All samples MUST have at least 5 average TPM.
ddF <- lowExpression_filter(ddT, groupsFactor_Trasnl, thres = 5, samples = 6)

lightHavgT <- rowMeans(ddF[,1:3])
lightLavgT <- rowMeans(ddF[,4:6])
heavyHavgT <- rowMeans(ddF[,7:9])
heavyLavgT <- rowMeans(ddF[,10:12])
totalHavgT <- rowMeans(ddF[,13:15])
totalLavgT <- rowMeans(ddF[,16:18])
tpmAvgAllF <- data.frame("LightL" = lightLavgT, "LightH" = lightHavgT, "HeavyL" = heavyLavgT, "HeavyH" = heavyHavgT, "TotalL" = totalLavgT, "TotalH" = totalHavgT)
tpmTrRatAllF <- data.frame("LightL" = lightLavgT/totalLavgT, "LightH" = lightHavgT/totalHavgT, "HeavyL" = heavyLavgT/totalLavgT, "HeavyH" = heavyHavgT/totalHavgT)

# Generate the Translation Ratio Data frame.
tpmTrRatF <- data.frame("HighGlu" = (tpmTrRatAllF$LightH + tpmTrRatAllF$HeavyH)/2, "LowGlu" = (tpmTrRatAllF$LightL + tpmTrRatAllF$HeavyL)/2, row.names = rownames(tpmTrRatAllF))

# Calculate the FC difference of translation ratio (normalised).
diffTranslRatioFC <- data.frame(transRatioFC = (tpmTrRatF$HighGlu - tpmTrRatF$LowGlu)/tpmTrRatF$LowGlu, row.names = rownames(tpmTrRatF))
diffTranslRatioFC_sorted <- diffTranslRatioFC[order(diffTranslRatioFC$transRatioFC, decreasing = TRUE), , drop = FALSE]

# Select the most (UP), the least (DOWN) and the non changing genes (Control).
diffTranslRatioFC_up <- rownames(diffTranslRatioFC_sorted[diffTranslRatioFC_sorted$transRatioFC >= 0.5, , drop = FALSE]) # 147
diffTranslRatioFC_down <- rownames(diffTranslRatioFC_sorted[diffTranslRatioFC_sorted$transRatioFC <= -0.25, , drop = FALSE]) # 154
diffTranslRatioFC_control <- rownames(diffTranslRatioFC_sorted[diffTranslRatioFC_sorted$transRatioFC >= -0.01 & diffTranslRatioFC_sorted$transRatioFC <= 0.01, , drop = FALSE]) # 326

# Keep all thse together in a file so that we will calculate the features.
write(c(diffTranslRatioFC_down, diffTranslRatioFC_up, diffTranslRatioFC_control), file = "diffTranslRatioFC_ALL.txt")

## calculate rna features with the rna_feat_ext suit. (627 genes were retrieved from 627) (the previous analysis with looser expression filter has retrieved 617 out of 717 genes)

# Load the features table.
featTranslRatioFC <- read.table("trans_ratio/diffTranslRatioFC_ALL.tab", sep = ";", header = TRUE)
# Keep only the numerical values.
featTranslRatioFC <- featTranslRatioFC[,1:14] # The first 14 columns contain the numerical values.

# Add the column for the translation ratio different. UP, Down, Control.
featTranslRatioFC["diffTranslRat"] <- NaN

for (i in 1:nrow(featTranslRatioFC)) {
  if (featTranslRatioFC[i, "ensembl_gene_id"] %in% diffTranslRatioFC_up) {
    featTranslRatioFC[i, "diffTranslRat"] = "UP"
  }
  else if (featTranslRatioFC[i, "ensembl_gene_id"] %in% diffTranslRatioFC_down) {
    featTranslRatioFC[i, "diffTranslRat"] = "DOWN"
  }
  else if (featTranslRatioFC[i, "ensembl_gene_id"] %in% diffTranslRatioFC_control) {
    featTranslRatioFC[i, "diffTranslRat"] = "Control"
  }
}

# Make the plots.
# Plot the Coding length boxplots
ptCD <- ggplot(featTranslRatioFC[featTranslRatioFC$coding_len <= 5000,], aes(x = diffTranslRat, y = coding_len, fill = diffTranslRat)) +
  theme_pubr() + scale_fill_brewer(palette="Dark2") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.15, height = 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  ylab("Coding Length") + xlab("Transl. Groups")
#  ggtitle("Coding length distribution of the 3 transl. difference genes sets") +
#  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill")

ptCDa <- ggbetweenstats( #Group clusters
  ggtheme = theme_pubr(),
  data = featTranslRatioFC[featTranslRatioFC$coding_len <= 6000,],
  ylab = "Coding Length", xlab = "H/L glucose translation ratio difference",
  x = diffTranslRat, y = coding_len,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  #title = "Coding length distribution of the 3 groups of translation differences.",
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE, mean.ci = TRUE) + ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) # Set FALSE for the manuscript figure.



# Plot GC content boxplots
ptGC <-ggplot(featTranslRatioFC, aes(x = diffTranslRat, y = GC, fill = diffTranslRat)) +
  theme_pubr() + scale_fill_brewer(palette="Dark2") +
  coord_cartesian(ylim = c(30, 75)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.15, height = 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  ylab("GC") + xlab("Transl. Groups")
#  ggtitle("Coding length distribution of the 3 transl. difference genes sets") +
#  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill")

ptGCa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featTranslRatioFC,
  x = diffTranslRat, y = GC,
  ylab = "GC", xlab = "H/L glucose translation ratio difference",
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 2, mean.point.args = list(size = 1.5, color = "darkred"),
  #title = "GC content of the 3 groups of translation differences.",
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE, mean.ci = TRUE) + ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) # Set FALSE for the manuscript figure.

# Plot 5'UTR length boxplots
pt5plen <- ggplot(featTranslRatioFC, aes(x = diffTranslRat, y = X5pUTR_len, fill = diffTranslRat)) +
  theme_pubr() + scale_fill_brewer(palette="Dark2") +
  coord_cartesian(ylim = c(0, 1000)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.15, height = 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  ylab("5'UTR length") + xlab("Transl. Groups")
#  ggtitle("Coding length distribution of the 3 transl. difference genes sets") +
#  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill")

pt5plena <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featTranslRatioFC[featTranslRatioFC$X5pUTR_len <= 1000 & featTranslRatioFC$X5pUTR_len != 0,],
  x = diffTranslRat, y = X5pUTR_len,
  ylab = "5'UTR length", xlab = "H/L glucose translation ratio difference",
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  #title = "5'UTR length of the 3 groups of translation differences.",
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE, mean.ci = TRUE) + ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) # Set FALSE for the manuscript figure.


# Plot 5'UTR GC boxplots
pt5pGC <- ggplot(featTranslRatioFC[featTranslRatioFC$X5pUTR_GC != 0,], aes(x = diffTranslRat, y = X5pUTR_GC, fill = diffTranslRat)) +
  theme_pubr() + scale_fill_brewer(palette="Dark2") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.15, height = 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  ylab("5'UTR GC") + xlab("Transl. Groups")
#  ggtitle("Coding length distribution of the 3 transl. difference genes sets") +
#  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill")

pt5pGCa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featTranslRatioFC[featTranslRatioFC$X5pUTR_GC != 0,],
  ylab = "5'UTR GC", xlab = "H/L glucose translation ratio difference",
  x = diffTranslRat, y = X5pUTR_GC,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 2, mean.point.args = list(size = 1.5, color = "darkred"),
  #title = "5'UTR GC of the 3 groups of translation differences.",
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE, mean.ci = TRUE) + ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) # Set FALSE for the manuscript figure.


# Plot 5'UTR MFE boxplots
pt5pMFE <- ggplot(featTranslRatioFC[featTranslRatioFC$X5pUTR_MFE != 0,], aes(x = diffTranslRat, y = X5pUTR_MFE, fill = diffTranslRat)) +
  theme_pubr() + scale_fill_brewer(palette="Dark2") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.15, height = 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  ylab("5'UTR MFE") + xlab("Transl. Groups")
#  ggtitle("Coding length distribution of the 3 transl. difference genes sets") +
#  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill")

pt5pMFEa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featTranslRatioFC[featTranslRatioFC$X5pUTR_MFE != 0,],
  ylab = "5'UTR MFE", xlab = "H/L glucose translation ratio difference",
  x = diffTranslRat, y = X5pUTR_MFE,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  #title = "5'UTR MFE of the 3 groups of translation differences.",
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE, mean.ci = TRUE) + ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) # Set FALSE for the manuscript figure.


# Plot 5'UTR MFE per BP boxplots
pt5pMfeBP <- ggplot(featTranslRatioFC[featTranslRatioFC$X5pUTR_MfeBP != 0,], aes(x = diffTranslRat, y = X5pUTR_MfeBP, fill = diffTranslRat)) +
  theme_pubr() + scale_fill_brewer(palette="Dark2") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.15, height = 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  ylab("5'UTR MFE/Bp") + xlab("Transl. Groups")
#  ggtitle("Coding length distribution of the 3 transl. difference genes sets") +
#  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill")

pt5pMfeBPa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featTranslRatioFC[featTranslRatioFC$X5pUTR_MfeBP != 0,],
  ylab = "5'UTR MFE per BP", xlab = "H/L glucose translation ratio difference",
  x = diffTranslRat, y = X5pUTR_MfeBP,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 2, mean.point.args = list(size = 1.5, color = "darkred"),
  #title = "5'UTR MFE per BP of the 3 groups of translation differences.",
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE, mean.ci = TRUE) + ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) # Set FALSE for the manuscript figure.


# Plot 3'UTR length boxplots
pt3pLEN <- ggplot(featTranslRatioFC[featTranslRatioFC$X3pUTR_len <= 5000 & featTranslRatioFC$X3pUTR_len != 0,], aes(x = diffTranslRat, y = X3pUTR_len, fill = diffTranslRat)) +
  theme_pubr() + scale_fill_brewer(palette="Dark2") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.15, height = 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  ylab("3'UTR length") + xlab("Transl. Groups")
#  ggtitle("3'UTR length distributions of the 3 translatin difference genes sets") +
#  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill")

pt3pLENa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featTranslRatioFC[featTranslRatioFC$X3pUTR_len <= 5000 & featTranslRatioFC$X3pUTR_len != 0,],
  x = diffTranslRat, y = X3pUTR_len,
  ylab = "3'UTR length", xlab = "H/L glucose translation ratio difference",
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  #title = "3'UTR length of the 3 groups of translation differences.",
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE, mean.ci = TRUE) + ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) # Set FALSE for the manuscript figure.


# Plot 3'UTR GC boxplots
pt3pGC <- ggplot(featTranslRatioFC[featTranslRatioFC$X3pUTR_GC != 0,], aes(x = diffTranslRat, y = X3pUTR_GC, fill = diffTranslRat)) +
  theme_pubr() + scale_fill_brewer(palette="Dark2") +
  coord_cartesian(ylim = c(15, 80)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.15, height = 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  ylab("3'UTR GC") + xlab("Transl. Groups")
#  ggtitle("3'UTR length distributions of the 3 translatin difference genes sets") +
#  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill")

pt3pGCa <- ggbetweenstats( # Group clusters
  ggtheme = theme_pubr(),
  data = featTranslRatioFC[featTranslRatioFC$X3pUTR_GC != 0,],
  ylab = "3'UTR GC", xlab = "H/L glucose translation ratio difference",
  x = diffTranslRat, y = X3pUTR_GC,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 2, mean.point.args = list(size = 1.5, color = "darkred"),
  #title = "3'UTR GC of the 3 groups of translation differences.",
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE, mean.ci = TRUE) + ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) # Set FALSE for the manuscript figure.


# Plot 3'UTR MFE boxplots
pt3pMFE <- ggplot(featTranslRatioFC[featTranslRatioFC$X3pUTR_MFE >= -2500 & featTranslRatioFC$X3pUTR_MFE != 0,], aes(x = diffTranslRat, y = X3pUTR_MFE, fill = diffTranslRat)) +
  theme_pubr() + scale_fill_brewer(palette="Dark2") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.15, height = 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  ylab("3'UTR MFE") + xlab("Transl. Groups")
#  ggtitle("3'UTR length distributions of the 3 translatin difference genes sets") +
  #stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill")

pt3pMFEa <- ggbetweenstats(
  ggtheme = theme_pubr(),
  data = featTranslRatioFC[featTranslRatioFC$X3pUTR_MFE >= -2500 & featTranslRatioFC$X3pUTR_MFE != 0,],
  ylab = "3'UTR MFE", xlab = "H/L glucose translation ratio difference",
  x = diffTranslRat,y = X3pUTR_MFE,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  #title = "3'UTR MFE of the 3 groups of translation differences.",
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE, mean.ci = TRUE) + ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) # Set FALSE for the manuscript figure.


# Plot 3'UTR MFE/BP boxplots
pt3pMfeBP <- ggplot(featTranslRatioFC[featTranslRatioFC$X3pUTR_MfeBP != 0,], aes(x = diffTranslRat, y = X3pUTR_MfeBP, fill = diffTranslRat)) +
  #coord_cartesian(ylim = c(0, -0.5)) +
  theme_pubr() + scale_fill_brewer(palette="Dark2") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.15, height = 0)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  ylab("3'UTR MFE/Bp") + xlab("Transl. Groups")
#  ggtitle("3'UTR length distributions of the 3 translatin difference genes sets") +
#  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill")

pt3pMfeBPa <- ggbetweenstats(
  ggtheme = theme_pubr(),
  data = featTranslRatioFC[featTranslRatioFC$X3pUTR_MfeBP != 0,],
  ylab = "3'UTR MFE per BP", xlab = "H/L glucose translation ratio difference",
  x = diffTranslRat, y = X3pUTR_MfeBP,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 2, mean.point.args = list(size = 1.5, color = "darkred"),
  #title = "3'UTR MFE per BP of the 3 groups of translation differences.",
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE, mean.ci = TRUE) + ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) # Set FALSE for the manuscript figure.


# Plot TOP local score boxplots
ptTOP <- ggplot(featTranslRatioFC, aes(x = diffTranslRat, y = TOP_localScore, fill = diffTranslRat)) +
  #coord_cartesian(ylim = c(0, -0.6)) +
  theme_pubr() + scale_fill_brewer(palette="Dark2") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.15, height = 0.05)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  ylab("TOP local score") + xlab("Transl. Groups")
#  ggtitle("3'UTR length distributions of the 3 translatin difference genes sets") +
#  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill")

ptTOPa <- ggbetweenstats(
  ggtheme = theme_pubr(),
  data = featTranslRatioFC,
  ylab = "TOP local score", xlab = "H/L glucose translation ratio difference",
  x = diffTranslRat, y = TOP_localScore,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 0, mean.point.args = list(size = 1.5, color = "darkred"),
  #title = "TOP local score of the 3 groups of translation differences.",
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE, mean.ci = TRUE) + ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) # Set FALSE for the manuscript figure.

# Plot CAI boxplots
ptCAI <- ggplot(featTranslRatioFC, aes(x = diffTranslRat, y = CAI, fill = diffTranslRat, group = diffTranslRat)) +
  #coord_cartesian(ylim = c(0, -0.6)) +
  theme_pubr() + scale_fill_brewer(palette="Dark2") +
  geom_boxplot(varwidth = TRUE, alpha = 0.5, notch = TRUE, outlier.shape = NA) +
  geom_jitter(alpha=0.4,  shape = 16, color = "grey2", position = position_jitter(width = 0.15, height = 0.05)) +
  geom_boxplot(varwidth = TRUE, alpha = 0.4, notch = TRUE, outlier.shape = NA) +
  theme(legend.position = "none", axis.title.x=element_blank()) +
  ylab("CAI") + xlab("Transl. Groups")
#  ggtitle("3'UTR length distributions of the 3 translatin difference genes sets") +
#  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5, position = "fill")

ptCAIa <- ggbetweenstats(
  ggtheme = theme_pubr(),
  data = featTranslRatioFC,
  ylab = "C Adaptation Index", xlab = "H/L glucose translation ratio difference",
  x = diffTranslRat, y = CAI,
  notch = TRUE, point.jitter.width = 1,
  type = "np", conf.level = 0.95, var.equal = FALSE,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.2), alpha = 0.5, size = 4, stroke = 0),
  k = 3, mean.point.args = list(size = 1.5, color = "darkred"),
  #title = "CAI of the 3 groups of translation differences.",
  pairwise.comparisons = TRUE, results.subtitle = TRUE, sample.size.label = TRUE, mean.plotting = TRUE, mean.ci = TRUE) + ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) # Set FALSE for the manuscript figure.

translEffBoxes <- ggarrange(ptCD, ptGC, pt5plen, pt5pGC, pt5pMFE, pt5pMfeBP,
                            pt3pLEN, pt3pGC, pt3pMFE, pt3pMfeBP, ptTOP, ptCAI,
                            labels = c("A", "B", "C", "D", "E" ,"F", "G", "H",
                                       "I", "J", "K", "L"),
                            ncol = 2, nrow = 6)

translEffBoxesAnal <- ggarrange(ptCDa, ptGCa, pt5plena, pt5pGCa, pt5pMFEa,
                                pt5pMfeBPa, pt3pLENa, pt3pGCa, pt3pMFEa,
                                pt3pMfeBPa, ptTOPa, ptCAIa,
                                labels = c("A", "B", "C", "D", "E" ,"F", "G",
                                           "H", "I", "J", "K", "L"),
                                ncol = 2, nrow = 6)


## TOP RNA analysis ---------------------------------------
# Known TOPrna from the union of 3 publications.
topRNAs <- read.table("publishedTOPrnas.csv", header=TRUE, sep = "\t")
topRNAsList <- list("Known" = as.vector(topRNAs$union.of.3[topRNAs$union.of.3 != ""]), "Cluster1" = as.vector(topRNAs$DEGs_Cluster.1[topRNAs$DEGs_Cluster.1 != ""]), "Cluster6" = as.vector(topRNAs$DEGs_Cluster.6[topRNAs$DEGs_Cluster.6 != ""]))
plot(euler(topRNAsList, shape = "ellipse"), quantities = TRUE)


## Further TOP RNA analysis --------------------------------
# Generate the data frame for the analysis of cluster1 and cluster6
#!!! WE SELECT top local score and the coding lenght as the first features, however we will explore more regressions (all vs. all).
featDF_1_6 <- featDF[(featDF$Cluster == 1 | featDF$Cluster == 6), ]
# generate the "Known" vector.
featDF_1_6$Known <- FALSE
for (r in 1:nrow(featDF_1_6)) {
  row = featDF_1_6[r,,]
  if (as.character(row$gene_name) %in% topRNAsList$Known) {
    featDF_1_6[r,]$Known <- TRUE
  }
}

# Generate a data frame by removing some outliers (i.e.e long and zero UTRs).
featDFtrim <- subset(featDF, featDF$coding_len <= 5000 & featDF$len_5pUTR <= 1000 & featDF$len_3pUTR <= 5000 & featDF$len_5pUTR != 0 & featDF$len_3pUTR != 0 & featDF$MfeBP_5pUTR != 0 & featDF$MfeBP_3pUTR != 0)

ggplot(subset(featDF_1_6, featDF_1_6$Cluster == 1), aes(x = TOP_localScore, y = log2(coding_len), color = Known)) + geom_point(size = 3, alpha = 0.75) + geom_rug() + ggtitle("Cluster 1")
ggplot(subset(featDF_1_6, featDF_1_6$Cluster == 6), aes(x = TOP_localScore, y = log2(coding_len), color = Known)) + geom_point(size = 3, alpha = 0.75) + geom_rug() + ggtitle("Cluster 6")

# All vs all correlation scatterplots for clusters.
ggpairs(subset(featDF_1_6[,3:14], featDF_1_6$Cluster == 1))
ggpairs(subset(featDF_1_6[,3:14], featDF_1_6$Cluster == 6))  # Simple ones.
ggpairs(featDF_1_6[rownames(featDFtrim),], columns = c(3, 4, 5, 6, 8, 9, 10, 12, 13 ,14), lower = list(continuous = "smooth"), diag =  list(continuous = "densityDiag"), upper = list(continuous = "density"), ggplot2::aes(colour=Cluster))  # Combined, trimmed etc.

# allVSall for the whole set as translation 1+6, 2+3, 4+5
ggpairs(featDFtrim, columns = c(3, 4, 5, 6, 8, 9, 10, 12, 13 ,14), lower = list(continuous = "smooth"), diag =  list(continuous = "densityDiag"), upper = list(continuous = "cor"), ggplot2::aes(colour=Translation, alpha = 0.5))
ggpairs(featDFtrim, columns = c(3, 4, 5, 6, 8, 9, 10, 12, 13 ,14), lower = list(continuous = "smooth"), diag =  list(continuous = "densityDiag"), upper = list(continuous = "density"), ggplot2::aes(colour=Translation, alpha = 0.5))


## UTRDB analysis -----------------------------------------
# Read the table from the UTRDB website analysis (the data file is preprocessed!!!!)
utr_5_table <- read.table("rna_feat/201904/utrScan_5utr_results.txt", sep = ":", header = TRUE, strip.white = TRUE)
utr_3_table <- read.table("rna_feat/201904/utrScan_3utr_results.txt", sep = ":", header = TRUE, strip.white = TRUE)

# Keep only the features column
utr_5_slice <- as.data.frame(utr_5_table[, c("trascript_ID", "feature")])
utr_3_slice <- as.data.frame(utr_3_table[, c("trascript_ID", "feature")])

# Do the count of each feature per transcript.
utr5 <- as.data.frame(utr_5_slice %>% dplyr::count(trascript_ID, feature, sort = TRUE, .drop = FALSE))
utr3 <- as.data.frame(utr_3_slice %>% dplyr::count(trascript_ID, feature, sort = TRUE, .drop = FALSE))

# Spread the data.
utr5_matrix_features <- utr5 %>%  spread(key = feature, value = n, fill = 0)
utr3_matrix_features <- utr3 %>%  spread(key = feature, value = n, fill = 0)

# Sort out the rownames
rownames(utr5_matrix_features) <- utr5_matrix_features$trascript_ID
utr5_matrix_features$trascript_ID <- NULL
rownames(utr3_matrix_features) <- utr3_matrix_features$trascript_ID
utr3_matrix_features$trascript_ID <- NULL

# Preprocess add an extra column from the clustering
utr5_matrixPlot <- utr5_matrix_features
utr5_matrixPlot$Clust <- 0
for (t in rownames(utr5_matrixPlot)){
  tt <- subset(featDF, featDF$X == t)
  utr5_matrixPlot[t,]$Clust <- tt$Cluster
}
dist5pUTR <- dist(utr5_matrixPlot[,1:3], method = "manhattan")  # Not used!
clust5pUTR <- hclust(dist5pUTR, method = "ward.D")  # Not used!
# Add the common gene names instead of ENSEMBLE IDs.
featDFi <- featDF
utr5_matrixPlot$gene_name <- featDFi[rownames(utr5_matrixPlot),, ]$gene_name
rownames(utr5_matrixPlot) <- utr5_matrixPlot$gene_name
utr5_matrixPlot$gene_name <- NULL
utr5_matrixPlot$PAS <- NULL

# For the 5'UTRs we have a range between 0 and 13 that's why we use 14 colours.
#colours5p <-  c("white", "#C6DBEF", "#4292C6", replicate(4, "#2171B5"), replicate(7, "#08306B"))
colfunc <- colorRampPalette(c("white", "dodgerblue4"))
colours5p <- colfunc(14)
colours5p <- c("#FFFFFF", "#A3BAD2", "#6B92B7", rep("#34699C", 4), rep("#104E8B", 7))
# 5'UTRs plotting clustering
hclust.local <- function(x) hclust(x, method="average")
par(mar = c(0,0,0,0))
heatmap.2(as.matrix(utr5_matrixPlot[,!names(utr5_matrixPlot) %in% "Clust"]), scale = "none", col = colours5p, RowSideColors = brewer.pal(n = 6, name = "Dark2")[as.factor(as.character(utr5_matrixPlot$Clust))], margins = c(9, 8), cexRow = 1.2, cexCol = 2.5, hclustfun = hclust.local, trace = "none", lhei = c(0.00001, 4), lwid = c(1, 5), key = FALSE, srtCol = 40)
legend("topleft", legend=c(1:6), fill= brewer.pal(6, "Dark2"), title = "Clusters", border = brewer.pal(6, "Dark2"), bty = "n", cex = 2, y.intersp = 1.5)
# Keep the informative part of the clustering
utr5_matrixPlotI <- utr5_matrixPlot[rowSums(utr5_matrixPlot[,c(10,8,3)] > 0) != 0,][, c(3,8,10,11)]
hclust.local <- function(x) hclust(x, method="ward.D2")
par(mar = c(0,0,0,0))
heatmap.2(as.matrix(utr5_matrixPlotI[,!names(utr5_matrixPlotI) %in% "Clust"]), scale = "none", col = colours5p, trace = "none", RowSideColors = brewer.pal(n = 6, name = "Dark2")[as.factor(as.character(utr5_matrixPlotI$Clust))], hclustfun = hclust.local, key = FALSE, srtCol = 45, margins = c(8, 7), lhei = c(0.00001, 4), lwid = c(1, 5), cexRow = 1.3, cexCol = 4)
legend("topleft", legend=c(1:6), fill = brewer.pal(6, "Dark2"), border = brewer.pal(6, "Dark2"), cex = 2, bty = "n", y.intersp = 1.5, title = "Cluster")

# Heatmap based on the model MClust ordering.
utr5_matrixPlotC <- utr5_matrixPlot[order(utr5_matrixPlot$Clust),]
heatmap.2(as.matrix(utr5_matrixPlotC[,names(utr5_matrixPlotC) %in% c("IRES", "TOP", "uORF")]), Colv = FALSE, Rowv = FALSE, dendrogram = "none", scale = "none", col = colours5p, RowSideColors = brewer.pal(n = 6, name = "Dark2")[as.factor(as.character(utr5_matrixPlotC$Clust))], margins = c(9, 8), cexRow = 1.2, cexCol = 2.5, trace = "none", lhei = c(0.00001, 4), lwid = c(1, 5), key = FALSE, srtCol = 45)  # NOT informative
legend("topleft", legend=c(1:6), fill = brewer.pal(6, "Dark2"), border = brewer.pal(6, "Dark2"), cex = 2, bty = "n", y.intersp = 1.5, title = "Cluster")


# Select the clusters of  mRNAs with a uORF
uorfDF <- data.frame("Cluster" = utr5_matrixPlot[rownames(subset(utr5_matrixPlot, utr5_matrixPlot$uORF>0)),]$Clust, row.names = rownames(subset(utr5_matrixPlot, utr5_matrixPlot$uORF>0)))
uorfDF <- uorfDF[order(uorfDF$Cluster, decreasing = TRUE),, drop = FALSE]


# 3'UTRs scan plotting clustering
# Preprocess add an extra column from the clustering
utr3_matrixPlot <- utr3_matrix_features
utr3_matrixPlot$Clust <- 0
for (t in rownames(utr3_matrixPlot)){
  tt <- subset(featDF, featDF$X == t)
  utr3_matrixPlot[t,]$Clust <- tt$Cluster
}
dist3pUTR <- dist(utr3_matrixPlot[,1:3], method = "manhattan")  # Not used
clust3pUTR <- hclust(dist3pUTR, method = "ward.D")  # Not used
# Add the common gene names instead of ENSEMBLE IDs.
featDFi <- featDF
rownames(featDFi) <- featDFi$X
featDFi$X <- NULL
utr3_matrixPlot$gene_name <- featDFi[rownames(utr3_matrixPlot), ]$gene_name
rownames(utr3_matrixPlot) <- utr3_matrixPlot$gene_name
utr3_matrixPlot$gene_name <- NULL

# For the 3'UTRs we have a range between 0 and 29 that's why we use 30 colours.
colfunc <- colorRampPalette(c("white", "dodgerblue4"))
colours3p <- colfunc(30)
colours3p <- c("#FFFFFF", "#A4BBD3", rep("#7B9DBF", 2), rep("#628BB3", 4), rep("#4978A7", 8), rep("#286096", 10), rep("#104E8B", 4))
# 3'UTRs plotting clustering
hclust.local <- function(x) hclust(x, method="ward.D")
par(mar = c(0,0,0,0))
heatmap.2(as.matrix(utr3_matrixPlot[,!names(utr3_matrixPlot) %in% "Clust"]), scale = "none", col = colours3p, RowSideColors = brewer.pal(n = 6, name = "Dark2")[as.factor(as.character(utr3_matrixPlot$Clust))], margins = c(10, 8), cexRow = 1.2, cexCol = 2.5, hclustfun = hclust.local, trace = "none", lhei = c(0.00001, 4), lwid = c(1, 5), key = FALSE, srtCol = 45)
legend("topleft", legend=c(1:6), fill= brewer.pal(6, "Dark2"), title = "Clusters", border = brewer.pal(6, "Dark2"), bty = "n", cex = 2, y.intersp = 1.5)
# one more time...
heatmap.2(as.matrix(utr3_matrixPlot[,!names(utr3_matrixPlot) %in% c("Clust", "PAS", "BRE", "GLUT1", "INS_SCE", "TGE", "IRE", "SECIS1", "SECIS2", "15-LOX-DICE")]), scale = "none", col = colours3p, RowSideColors = brewer.pal(n = 6, name = "Dark2")[as.factor(as.character(utr3_matrixPlot$Clust))], margins = c(9, 8), cexRow = 1.2, cexCol = 2.5, hclustfun = hclust.local, trace = "none", lhei = c(0.00001, 4), lwid = c(1, 5), key = FALSE, srtCol = 45)
legend("topleft", legend=c(1:6), fill= brewer.pal(6, "Dark2"), title = "Clusters", border = brewer.pal(6, "Dark2"), bty = "n", cex = 2, y.intersp = 1.5)
# Heatmap based on the model MClust ordering.
utr3_matrixPlotC <- utr3_matrixPlot[order(utr3_matrixPlot$Clust),]
heatmap.2(as.matrix(utr3_matrixPlotC[,!names(utr3_matrixPlotC) %in% c("Clust", "PAS", "BRE", "GLUT1", "INS_SCE", "TGE", "IRE", "SECIS1", "SECIS2", "15-LOX-DICE")]), Colv = FALSE, Rowv = FALSE, dendrogram = "none", scale = "none", col = colours3p, RowSideColors = brewer.pal(n = 6, name = "Dark2")[as.factor(as.character(utr3_matrixPlotC$Clust))], margins = c(9, 8), cexRow = 1.2, cexCol = 2.5, trace = "none", lhei = c(0.00001, 4), lwid = c(1, 5), key = FALSE, srtCol = 45)  # NOT informative

# Keep the informative part of the clustering  NOT so informative for the 3'UTR
utr3_matrixPlotI <- utr3_matrixPlot[rowSums(utr3_matrixPlot[,c("K-BOX", "SXL_BS", "IRES", "GY-BOX", "UNR-bs", "BRD-BOX", "CPE")] > 0) != 0,][,c("K-BOX", "SXL_BS", "IRES", "GY-BOX", "UNR-bs", "BRD-BOX", "CPE", "Clust")]
hclust.local <- function(x) hclust(x, method="average")
par(mar = c(0,0,0,0))
heatmap.2(as.matrix(utr3_matrixPlotI[,!names(utr3_matrixPlotI) %in% "Clust"]), scale = "none", col = colours3p, trace = "none", RowSideColors = brewer.pal(n = 6, name = "Dark2")[as.factor(as.character(utr3_matrixPlotI$Clust))], hclustfun = hclust.local, key = FALSE, srtCol = 45, margins = c(8, 7), lhei = c(0.00001, 4), lwid = c(1, 5), cexRow = 1.25, cexCol = 2.5)
legend("topleft", legend=c(1:6), fill = brewer.pal(6, "Dark2"), border = brewer.pal(6, "Dark2"), cex = 2, bty = "n", y.intersp = 1.5, title = "Cluster")


utr3_features_final <- utr3_matrix_features[-c(20,1,5,7,9,10,13,17,18,3,15,14)]
row_sub = apply(utr3_features_final, 1, function(row) any(row != 0 ))
utr3_features_final <- utr3_features_final[row_sub,]
colours3p <-  c("white", "#C6DBEF", "#4292C6", "#2171B5", "#08306B")
heatmap(as.matrix(utr3_matrix_features), scale = "none", col = colours3p)
heatmap.2(as.matrix(utr3_features_final[,c(14,11,4,8,16,6,2)]), scale = "none", col = colours3p, Colv = FALSE, dendrogram = "row", trace = "none", key = FALSE, margins = c(7, 7))


#G4 analysis
G4_genes <- read.table("rna_feat/201904/G4_5UTR_names.txt", header = F, sep = ";")
G4_unique <- unique(G4_genes$V1)
degs_g4_names <- intersect(G4_unique, degs)
G4_names_cluster <- clusterGeneIDs[degs_g4_names,, drop = FALSE]
hist(G4_names_cluster$cluster, nclass = 50)

freq_g4_cluster <- as.data.frame((table(G4_names_cluster$cluster)))
freq_cluster <- as.data.frame((table(clusterGeneIDs$cluster)))

freq_cluster_g4cluster <- cbind(freq_cluster, freq_g4_cluster$Freq)
colnames(freq_cluster_g4cluster) <- c("ClusterID","Size","G4_size")

freq_cluster_g4cluster <- freq_cluster_g4cluster %>% mutate(Percentage = (G4_size/Size)*100)
hist(freq_cluster_g4cluster$Percentage,nclass = 25)

translation <- c("Up","Inter","Inter", "Down", "Down", "Up")
freq_cluster_g4cluster$translation <- translation

ggplot(freq_cluster_g4cluster, aes(x = ClusterID, y = Percentage, fill = translation, group = ClusterID)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  theme_pubr() +
  scale_fill_brewer(palette="Dark2") +
  theme(text = element_text(size = 22)) +
  xlab("Cluster ID")
  #+ ggtitle("Percentrage of G-quadruplex genes of the 6 Mclust clusters")


# REVISIONS --------

## EndoBeta microarray ---------------------------
transEndoB <- read.table("data/GSE48101_full_matrix.txt", header = TRUE, sep = "\t", colClasses = c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))

groupsEndoB <- factor(c("GSM99", "GSM99", "GSM99", "H357cre", "H357cre", "H357cre", "H357p59", "H357p59", "H357p59", "SKPC", "SKPC", "SKPC"), levels = c("GSM99", "H357cre", "H357p59", "SKPC"))

groupsEndoBB <- factor(c("H357cre", "H357cre", "H357cre", "H357p59", "H357p59", "H357p59", "SKPC", "SKPC", "SKPC"), levels = c("H357cre", "H357p59", "SKPC"))


# A bit of filtering perhaps...
dd <- transEndoB
dd$probeID <- rownames(transEndoB)
dd <- melt(dd, id.vars = "probeID")
dd$value <- as.numeric(dd$value)
dd$value <- log2(dd$value)
# Just
ggplot(dd, aes(x = value)) + facet_wrap(~variable,scales = "free_x") + geom_histogram()

# PCA on the microarray values directly.
res.pca.EndoB = PCA(t(transEndoB), graph = FALSE)
res.pca.EndoBB = PCA(t(transEndoB[,4:12]), graph = FALSE)

fviz_pca_ind(res.pca.EndoB,
             fill.ind = groupsEndoB, col.var = "black", repel = TRUE,
             col.ind = groupsEndoB, # coloured by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "orange"),
             addEllipses = TRUE, # Ellipses de concentration
             legend.title = "Groups",
             title = "PCA on the microarrays of endoBeta cells")

fviz_pca_ind(res.pca.EndoBB,
             fill.ind = groupsEndoBB, col.var = "black", repel = TRUE,
             col.ind = groupsEndoBB, # coloured by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Ellipses de concentration
             legend.title = "Groups",
             title = "PCA on the microarrays of endoBeta cells")

# Filter and assignee probes to genes.
affyGenesProbes <- read.table("data/genesProbes_HGU133_Plus_2.csv", header = TRUE, sep = ",")
# Calculate MAD on the two H357 cell lines.
## perhaps we want to calculate the MAD between H357 and SKPC also....
transEndoB$MAD <- apply(transEndoB[,4:9], 1, mad)
# Sort by affy_probe names and add gene names.
affyGenesProbes <- affyGenesProbes[order(affyGenesProbes$affy_ProbeID),]
transEndoB$geneID <- affyGenesProbes$gene_Name

# Generate two intermediate matrices to select for probes with the higest MAD.
x <- transEndoB[order(transEndoB$MAD, dedecreasing = TRUE),]
y <- x[which(!duplicated(x$geneID)),]

# Put the H357 and SKPC lines in a new proper data.frame
transH357 <- data.frame(y[,4:12], row.names = y$geneID)
# Remove a funny line with gene name "---"
transH357 <- transH357[-c(2),]

# Remove lowly expresses genes.
h357_medians <- rowMedians(as.matrix(transH357))
# Plot row median intensities.
hist(h357_medians, 100, col = "cornsilk1", freq = FALSE,
     main = "Histogram of the median intensities",
     border = "antiquewhite4",
     xlab = "Median intensities")

# Differential expression
designMat_E <- model.matrix(~0+groupsEndoBB)
colnames(designMat_E) <- levels(groupsEndoBB)

contr.matrix_E <- makeContrasts(
  Cre_vs_P59 = H357cre - H357p59,
  Cre_vs_SKPC = H357cre - SKPC,
  P59_vs_SKPC = H357p59 - SKPC,
  levels = colnames(designMat_E))

# DiffExp Endo H357 SKPC
vmE <- voom(transH357, designMat_E, plot = TRUE)
fitE <- lmFit(vmE, designMat_E)
vfitE <- contrasts.fit(fitE, contrasts = contr.matrix_E)
efitE <- eBayes(vfitE, robust = TRUE)
efE <- decideTests(efitE, p.value = 0.01, lfc = 1)
summary(efE)
plotSA(efitE, main = "SA plot for probes from H357")
plotMD(efitE, column = 1, status = efE[,1], main = colnames(efitE)[1], ylim = c( -1.5, 1.2))
tfitE <- treat(vfitE, lfc = log2(1.5))
ttE <- decideTests(tfitE)
summary(ttE)
plotMD(tfitE, column = 1, status = ttE[,1], main = colnames(tfitE)[1], ylim = c( -1.5, 1.2))
plotMD(tfitE, column = 2, status = ttE[,2], main = colnames(tfitE)[2], ylim = c( -2, 2))
plotMD(tfitE, column = 3, status = ttE[,3], main = colnames(tfitE)[3], ylim = c( -2, 2))


# Show with mixomics that the two cell lines are similar....
