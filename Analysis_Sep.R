library(ggplot2)
library(DESeq2)
library(reshape)
library(VennDiagram)
library(pheatmap)
library(limma)
library(Glimma)
library(edgeR)
library("FactoMineR")
library("factoextra")
library(gplots)
library(RColorBrewer)

# Function to plot expression of individual gene.
geneBoxplotCond <- function(matrix, name, experiments, treatments, jit_width = 0.1, point_size = 2, ...){
  # Experiments are the different fractions.
  # Treatment is High or Low glucose.
  ge <- data.frame(t(matrix[name,]));
  ge$exp <- experiments;
  ge$treat <- treatments;
  colnames(ge)[1] <- "TPM";
  p <- ggplot(ge, aes(exp, TPM));
  p + geom_jitter(aes(color = treat), width = jit_width, size = point_size) + ggtitle(name);
}


tpmall <- read.table("polysomeProfile_TPM_proteinCoding.csv", header = TRUE)
tpmall <- as.data.frame(tpmall[, c(7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,1,2,3,4,5,6)])

cpmall <- as.data.frame(cpm(countsTableRaw))
cpmall <- as.data.frame(cpmall[, c(7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,1,2,3,4,5,6)])
head(cpmall)

countsTableRaw <- read.table("countsTOTALS_CodingGenes.tsv", header = TRUE)
countsTableRawH <- countsTableRaw[,19:24]
countsTableRawL <- countsTableRaw[,13:18]
countsTableRawM <- countsTableRaw[,7:12]
countsTableRawT <- countsTableRaw[,1:6]

groupsformatrix_L <- factor(c("LightH","LightH","LightH","LightL","LightL","LightL"),  
                                levels = c("LightH","LightL"))

groupsformatrix_H <- factor(c("HeavyH","HeavyH","HeavyH","HeavyL","HeavyL",
                                  "HeavyL"),  
                                levels = c("HeavyH","HeavyL"))

groupsformatrix_M <- factor(c("MonoH","MonoH","MonoH","MonoL","MonoL","MonoL"),  
                                levels = c("MonoH","MonoL"))

groupsformatrix_T <- factor(c("TotalH","TotalH","TotalH","TotalL","TotalL","TotalL"),  
                            levels = c("TotalH","TotalL"))

treasAll <- as.factor(c("High", "High", "High", "Low", "Low", "Low"))

cpmall <- cpm(countsTableRaw)

groupsall <- factor(c("Mono","Mono","Mono","Mono","Mono","Mono",
                      "Light","Light","Light","Light","Light","Light",
                      "Heavy","Heavy","Heavy","Heavy","Heavy","Heavy",
                      "Total","Total","Total","Total","Total","Total"),
                    levels = c("Mono", "Light", "Heavy", "Total"));

treasAll <- as.factor(c(rep(c("High", "High", "High", "Low", "Low", "Low"), 4)));


cpmtmpH <- cpm(countsTableRawH)
cpmtmpM <- cpm(countsTableRawM)
cpmtmpL <- cpm(countsTableRawL)
cpmtmpT <- cpm(countsTableRawT)

cpmtmpL <- cpmtmpL[apply(cpmtmpL>1, 1, all),]
cpmtmpH <- cpmtmpH[apply(cpmtmpH>1, 1, all),]
cpmtmpM <- cpmtmpM[apply(cpmtmpM>1, 1, all),]
cpmtmpT <- cpmtmpT[apply(cpmtmpT>1, 1, all),]

cpmallL <- lowExpression_filter(cpmtmpL, groupsformatrix_L, thres = 4, samples = 2)
cpmallH <- lowExpression_filter(cpmtmpH, groupsformatrix_H, thres = 4, samples = 2)
cpmallM <- lowExpression_filter(cpmtmpM, groupsformatrix_M, thres = 4, samples = 2)
cpmallT <- lowExpression_filter(cpmtmpT, groupsformatrix_T, thres = 4, samples = 2)

countsTable_H <- countsTableRawH[rownames(cpmallH),]
countsTable_L <- countsTableRawL[rownames(cpmallL),]
countsTable_M <- countsTableRawM[rownames(cpmallM),]
countsTable_T <- countsTableRawT[rownames(cpmallT),]

designMat_H <- model.matrix(~0+groupsformatrix_H)
colnames(designMat_H) <- levels(groupsformatrix_H)

designMat_M <- model.matrix(~0+groupsformatrix_M)
colnames(designMat_M) <- levels(groupsformatrix_M)

designMat_L <- model.matrix(~0+groupsformatrix_L)
colnames(designMat_L) <- levels(groupsformatrix_L)

designMat_T <- model.matrix(~0+groupsformatrix_T)
colnames(designMat_T) <- levels(groupsformatrix_T)

contr.matrix_H <- makeContrasts(
  HeavyPoly_HvsL = HeavyH - HeavyL,
  levels = colnames(designMat_H))

contr.matrix_M <- makeContrasts(
  Monosomes_HvsL = MonoH - MonoL,
  levels = colnames(designMat_M))

contr.matrix_L <- makeContrasts(
  LightPoly_HvsL= LightH - LightL,
  levels = colnames(designMat_L))

contr.matrix_T <- makeContrasts(
  Total_HvsL = TotalH - TotalL,
  levels = colnames(designMat_T))

vmL <- voom(countsTable_L, designMat_L, plot = TRUE)
fitL <- lmFit(vmL, designMat_L)
vfitL <- contrasts.fit(fitL, contrasts = contr.matrix_L)
efitL <- eBayes(vfitL, robust = TRUE) # Playing with the parameters of ebayes makes no difference.
efL <- decideTests(efitL, p.value = 0.05, lfc= 0.5)
summary(efL)
plotSA(efitL)
tfitL <- treat(vfitL, lfc= log2(1.1))
ttL <- decideTests(tfitL)
summary(ttL)

vmM <- voom(countsTable_M, designMat_M, plot = TRUE)
fitM <- lmFit(vmM, designMat_M)
vfitM <- contrasts.fit(fitM, contrasts = contr.matrix_M)
efitM <- eBayes(vfitM, robust = TRUE) # Playing with the parameters of ebayes makes no difference.
efM <- decideTests(efitM, p.value = 0.05, lfc= 0.5)
summary(efM)
plotSA(efitM)
tfitM <- treat(vfitM, lfc= log2(1.1))
ttM <- decideTests(tfitM)
summary(ttM)

vmH <- voom(countsTable_H, designMat_H, plot = TRUE)
fitH <- lmFit(vmH, designMat_H)
vfitH <- contrasts.fit(fitH, contrasts = contr.matrix_H)
efitH <- eBayes(vfitH, robust = TRUE) # Playing with the parameters of ebayes makes no difference.
efH <- decideTests(efitH,p.value = 0.05, lfc = 0.5)
summary(efH)
plotSA(efitH)
tfitH <- treat(vfitH, lfc= log2(1.1))
ttH <- decideTests(tfitH)
summary(ttH)

vmT <- voom(countsTable_T, designMat_T, plot = TRUE)
fitT <- lmFit(vmT, designMat_T)
vfitT <- contrasts.fit(fitT, contrasts = contr.matrix_T)
efitT <- eBayes(vfitT, robust = TRUE) # Playing with the parameters of ebayes makes no difference.
efT <- decideTests(efitT,p.value = 0.05, lfc= 0.5)
summary(efT)
plotSA(efitT)
tfitT <- treat(vfitT, lfc= log2(1.15))
ttT <- decideTests(tfitT)
summary(ttT)

plotMD(efitM, column = 1, status = efM[,1], main = colnames(efitM)[1],ylim = c( -1.5, 1.5))
plotMD(efitL, column = 1, status = efL[,1], main = colnames(efitL)[1],ylim = c( -1.5, 1.5))
plotMD(efitH, column = 1, status = efH[,1], main = colnames(efitH)[1],ylim = c( -1.5, 1.5))
plotMD(efitT, column = 1, main = colnames(efitT)[1],ylim = c( -1.5, 1.5))

PCA(t(vmM$E))

barplot(colSums(countsTable_M), col = c("lightgreen", "lightgreen","lightgreen","lightblue", "lightblue", "lightblue"), main = "Library size Raw counts")
barplot(colSums(vmM$E), col = c("lightgreen", "lightgreen","lightgreen","lightblue", "lightblue", "lightblue"), main = "Normalised counts (CPM voom method)",ylab = "CPM")


DEH <- topTable(efitH, coef=1, n=Inf, p.value = 0.05, lfc = 0.5)
DEL <- topTable(efitL, coef=1, n=Inf, p.value=0.05, lfc = 0.5)
DEM <- topTable(efitM, coef=1, n=Inf, p.value=0.05, lfc = 0.5)
DET <- topTable(efitT, coef=1, n=Inf, p.value=0.05, lfc = 0.5)

degs <- rbind(DEH,DEL)

write.csv(DEH, "DEG_H.csv")
write.csv(DEL, "DEG_L.csv")

vennDiagram(efH,efL, include = "up", circle.col = c("turquoise", "salmon"), main = "Venn with t-FIT")

 
NamesH <-row.names(DEH)
NamesL <-row.names(DEL)

NewDEG <- union(NamesL,NamesH)

NamesHL <- list("heavy" = NamesH, "light" = NamesL);
NamesHL <- calculate.overlap(NamesHL)



barplot(myGO$ngenes,names.arg =  myGO$categories, 
        col = my_palette, horiz = TRUE, las=1,
        main = "Functional Classification", xlab = "Number of genes")

bar <- ggplot(data = gsea, aes(x = factor(gsea$ï..Gene.Set.Name, levels = gsea$ï..Gene.Set.Name), y = gsea$Genes.in.Overlap..k., fill= gsea$FDR.q.value)) +
  geom_bar(stat = 'identity')+ 
  scale_fill_gradient2(low=LtoM(100), mid='snow3', 
                       high=MtoH(100), space='Lab')

print(bar)


write(NamesH, "Gene_DEH.txt")
write(NamesL, "Gene_DEL.txt")
venn.plot <- venn.diagram(NamesHL,
                          filename = "Venn_heavy-light.tiff",
                          scaled = TRUE,
                          ext.text = TRUE,
                          fill = 	"dodgerblue",
                          ext.line.lwd = 1,
                          ext.dist = -1,
                          ext.length = 5,
                          ext.pos = -3,
                          inverted = TRUE,
                          cat.dist = 0.05,
                          cex = 1,
                          cat.cex = 1,
                          rotation.degree = 0,
                          main = "Overlap Heavy Vs Light",
                          main.cex = 3)



source("https://bioconductor.org/biocLite.R")
biocLite("GOSemSim")
biocLite("org.Hs.eg.db")
biocLite("clusterProfiler")
biocLite("AnnotationDbi")
biocLite("ggplot2")
biocLite("DESeq2")
biocLite("reshape")
biocLite("VennDiagram")
biocLite("pheatmap")
biocLite("limma")
biocLite("Glimma")
biocLite("edgeR")
biocLite("FactoMineR")
biocLite("factoextra")
biocLite("gplots")
biocLite("RColorBrewer")
biocLite("ggdendro")
biocLite("GO.db")
biocLite("anota2seq")

install.packages("msigdbr")
install.packages("ReactomePA")

# Load packages.
library(magrittr)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ReactomePA)
library(enrichplot)
library(RColorBrewer)
library(gplots)



library(AnnotationDbi)
library(stringi)
library(digest)
library(htmlwidgets)
library(bit)
library(Rcpp)
library(DOSE)
library(topGO)
library(ALL)
library(clusterProfiler)
library(anota2seq)
library(GOSemSim)

### CLUSTER GO ANALYSIS
#IMPORT TABLE

degcluster <-read.csv("degs_16102018_Clusters.csv", header = TRUE, sep = "\t")

geneIDs1 <- subset(degcluster, degcluster$Cluster==1)
geneIDs2 <- subset(degcluster, degcluster$Cluster==2)
geneIDs3 <- subset(degcluster, degcluster$Cluster==3)
geneIDs4 <- subset(degcluster, degcluster$Cluster==4)
geneIDs5 <- subset(degcluster, degcluster$Cluster==5)
geneIDs6 <- subset(degcluster, degcluster$Cluster==6)
geneIDs7 <- subset(degcluster, degcluster$Cluster==7)

geneIDs1 <- as.character(unique(geneIDs1$ensembl_gene_id))
geneIDs2 <- as.character(unique(geneIDs2$ensembl_gene_id))
geneIDs3 <- as.character(unique(geneIDs3$ensembl_gene_id))
geneIDs4 <- as.character(unique(geneIDs4$ensembl_gene_id))
geneIDs5 <- as.character(unique(geneIDs5$ensembl_gene_id))
geneIDs6 <- as.character(unique(geneIDs6$ensembl_gene_id))
geneIDs7 <- as.character(unique(geneIDs7$ensembl_gene_id))

egoC1 <- clusterProfiler::enrichGO(gene = geneIDs1, OrgDb = OrgDb, ont = "BP", pAdjustMethod = "BH",  pvalueCutoff  = 0.05, qvalueCutoff  = 0.1,  readable = TRUE, keyType = "ENSEMBL")
egoC2 <- clusterProfiler::enrichGO(gene = geneIDs2, OrgDb = OrgDb, ont = "BP", pAdjustMethod = "BH",  pvalueCutoff  = 0.05, qvalueCutoff  = 0.1,  readable = TRUE, keyType = "ENSEMBL")
egoC3 <- clusterProfiler::enrichGO(gene = geneIDs3, OrgDb = OrgDb, ont = "BP", pAdjustMethod = "BH",  pvalueCutoff  = 0.05, qvalueCutoff  = 0.1,  readable = TRUE, keyType = "ENSEMBL")
egoC4 <- clusterProfiler::enrichGO(gene = geneIDs4, OrgDb = OrgDb, ont = "BP", pAdjustMethod = "BH",  pvalueCutoff  = 0.05, qvalueCutoff  = 0.1,  readable = TRUE, keyType = "ENSEMBL")
egoC5 <- clusterProfiler::enrichGO(gene = geneIDs5, OrgDb = OrgDb, ont = "BP", pAdjustMethod = "BH",  pvalueCutoff  = 0.05, qvalueCutoff  = 0.1,  readable = TRUE, keyType = "ENSEMBL")
egoC6 <- clusterProfiler::enrichGO(gene = geneIDs6, OrgDb = OrgDb, ont = "BP", pAdjustMethod = "BH",  pvalueCutoff  = 0.05, qvalueCutoff  = 0.1,  readable = TRUE, keyType = "ENSEMBL")
egoC7 <- clusterProfiler::enrichGO(gene = geneIDs7, OrgDb = OrgDb, ont = "BP", pAdjustMethod = "BH",  pvalueCutoff  = 0.05, qvalueCutoff  = 0.1,  readable = TRUE, keyType = "ENSEMBL")

egosimpleC1 <- simplify(egoC1, cutoff=0.5, by="p.adjust", select_fun=min, measure = "Wang")
egosimpleC2 <- simplify(egoC2, cutoff=0.5, by="p.adjust", select_fun=min, measure = "Wang")
egosimpleC3 <- simplify(egoC3, cutoff=0.5, by="p.adjust", select_fun=min, measure = "Wang")
egosimpleC4 <- simplify(egoC4, cutoff=0.5, by="p.adjust", select_fun=min, measure = "Wang")
egosimpleC5 <- simplify(egoC5, cutoff=0.5, by="p.adjust", select_fun=min, measure = "Wang")
egosimpleC6 <- simplify(egoC6, cutoff=0.5, by="p.adjust", select_fun=min, measure = "Wang")
egosimpleC7 <- simplify(egoC7, cutoff=0.5, by="p.adjust", select_fun=min, measure = "Wang")

# Further analysis!
##
# GO categorisation
OrgDb <- org.Hs.eg.db::org.Hs.eg.db
ggo <- clusterProfiler::groupGO(gene = NewDEG, OrgDb = OrgDb, ont = "BP", level = 5, keyType = "ENSEMBL")
barplot(ggo, drop = TRUE, showCategory = 30)


# GO enrichment
ego <- clusterProfiler::enrichGO(gene = rownames(DEheavyall),
                                 OrgDb = OrgDb, ont = "BP", pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05, qvalueCutoff  = 0.1, 
                                 readable = TRUE, keyType = "ENSEMBL")
barplot(egoUP, drop = TRUE, showCategory = 50,
        title = "GO Biological Pathways",
        font.size = 8
        )

emapplot(ego)
cnetplot(ego, categorySize = "pvalue", foldChange = DEH$logFC)

bp2 <- simplify(ego, cutoff=0.5, by="p.adjust", select_fun=min, measure = "Wang")


egoM <- clusterProfiler::enrichGO(gene = rownames(DEheavyall), OrgDb = OrgDb, ont = "BP", pAdjustMethod = "BH",  pvalueCutoff  = 0.05, qvalueCutoff  = 0.1,  readable = TRUE, keyType = "ENSEMBL")
barplot(egoM, drop = TRUE, showCategory = 50,
        title = "GO Molecular Function",
        font.size = 8)
goplot(ego)


# KEGG enrichment.
keggDEheavy <- scan("keggIDs.txt", what = "character")
kk <- clusterProfiler::enrichKEGG(gene = ids$ENTREZID, organism = "hsa", pvalueCutoff = 1)
barplot(kk, showCategory = 8)

# Pathview of the insuline secretion pathway.
pathview_func(DEheavyallC, logFCcolumn = "logFC", pathway.id = "04911", out.suffix = "insulineSecre")




## New analysis with the 438 genes of Manuel.
deManu <- scan("final_DEgenes.txt", what = "character")
DEfinal_Manuel_tpm <- tpmall[deManu,]
DEfinal_Manuel_tpm <- tpmall1[NewDEG,]

# Calculate the log Ratios over TPMs.
#Calculate means of groups
monoHall <- rowMeans(DEfinal_Manuel_tpm[,1:3])
monoLall <- rowMeans(DEfinal_Manuel_tpm[,4:6])
lightHall <- rowMeans(DEfinal_Manuel_tpm[,7:9])
lightLall <- rowMeans(DEfinal_Manuel_tpm[,10:12])
heavyHall <- rowMeans(DEfinal_Manuel_tpm[,13:15])
heavyLall <- rowMeans(DEfinal_Manuel_tpm[,16:18])

genesDEavgall <- data.frame("monoHall" = monoHall, "monoLall" = monoLall, "lightHall" = lightHall,"lightLall" = lightLall,"heavyHall"= heavyHall, "heavyLall" = heavyLall, row.names = names(heavyHall))

# Calculate logRatios over means
ratioHall <- log2(genesDEavgall$heavyHall/genesDEavgall$heavyLall)
ratioLall <- log2(genesDEavgall$lightHall/genesDEavgall$lightLall)
ratioMall <- log2(genesDEavgall$monoHall/genesDEavgall$monoLall)

logRatiosDEGtableManuel <- data.frame("Mono" = ratioMall, "Light" = ratioLall, "Heavy" = ratioHall,  row.names = rownames(genesDEavgall))

my_palette <- brewer.pal(n = 11, name = "RdYlGn")
# Clustering of the log ratios.

# Hierarchical clustering
heatmap.2(as.matrix(logRatiosDEGtableManuel), col = my_palette, dendrogram = "row", cexCol = 1.5, cexRow = 0.0001, key.title = NA, keysize = 0.9, main = "DEGs logRatio H/L hclust", ylab = "Genes")

## k-means clustering
# Perform k-means clustering k-means # 8 is empirically obtained from the 8 different possible behaviours of up-down
km <- kmeans(as.matrix(logRatiosDEGtableManuel), 8)

# Append id and cluster
dfcall <- cbind(logRatiosDEGtableManuel, id = seq(nrow(logRatiosDEGtableManuel)), cluster = km$cluster)

# Add idsort, the id number ordered by cluster
dfcall$idsort <- dfcall$id[order(dfcall$cluster)]
dfcall$idsort <- order(dfcall$idsort)

clusterCols <- as.character(sort(km$cluster))
#clusterCols <- str_replace_all(clusterCols, "1", "grey") #to change colors

heatmap(as.matrix(logRatiosDEGtableManuel)[order(km$cluster),], Rowv = NA, Colv = NA, scale = "none", labRow = NA, cexCol = 1.5, col = my_palette, RowSideColors = clusterCols, ylab = "Genes LogRatio", main = "k-means clustering of H/L LogRatios")





countsTableRaw <- read.table("countsTOTALS_CodingGenes.tsv", header = TRUE)
count_T <- countsTableRaw[, 1:6]
count_H <- countsTableRaw[, 19:24]

phenoVec <- as.factor(c("High", "High", "High", "Low", "Low", "Low"))

ads <- anota2seqDataSetFromMatrix(
  dataP = count_H,
  dataT = count_T,
  phenoVec = phenoVec,
  dataType = "RNAseq",
  filterZeroGenes = TRUE,
  normalize = TRUE,
  transformation = "TMM-log2",
  varCutOff = NULL)


ads<- anota2seqPerformQC(Anota2seqDataSet = ads,
                          generateSingleGenePlots = TRUE)

ads <- anota2seqResidOutlierTest(ads)

ads <- anota2seqAnalyze(Anota2seqDataSet = ads,
                        analysis = c("translation", "buffering"))

head(anota2seqGetOutput(
  ads, analysis = "translation",
  output = "full",
  selContrast = 1,
  getRVM = TRUE))

par(mfrow = c(1, 1))
anota2seqPlotPvalues(ads, selContrast = 1, plotToFile = FALSE)

ads <- anota2seqSelSigGenes(Anota2seqDataSet = ads,
                            selContrast = 1,
                            maxPAdj = 0.5)





ids <- bitr(NewDEG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
ids <- bitr(NewDEG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
