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


# Function to filter form lowly expressed genes.
lowExpression_filter <- function(e, groups, thres = 3, samples = 1, ...){
  # e is the raw counts data.frame (or better cpm matrix).
  # *groups* is a factor designated the goroupins (conditions etc.)
  ##it MUSt be of equal lenth to the columns of e.
  # thres is the threshold of the mean.
  # samples is the minimum number of sample groups that we want this threshold
  ##thres to be higher.
  filteredDF <- data.frame()
  rows <- vector();
  for (i in 1:nrow(e)) {
    row <- as.numeric(e[i,]);
    means <- aggregate(row~groups, FUN = mean);
    if (sum(means$row >= thres) >= samples) { # This condition checks for the grouping sums.
      filteredDF <- rbind(filteredDF, row);
      rows <- append(rows, rownames(e)[i])
    }
  }
  colnames(filteredDF) <- colnames(e);
  rownames(filteredDF) <- rows;
  return(filteredDF)
}


# Function to plot expression of individual gene.
geneBoxplotCond <- function(matrix, name, experiments, treatments, jit_width = 0.1, point_size = 2, ...){
  # Experiments are the different fractions.
  # Treatment is High or Low glucose.
  ge <- data.frame(t(matrix[name,]));
  ge$exp <- experiments;
  ge$treat <- treatments;
  colnames(ge)[1] <- "CPM";
  p <- ggplot(ge, aes(exp, CPM));
  p + geom_jitter(aes(color = treat), width = jit_width, size = point_size) + ggtitle(name);
}


# Get table of counts from STAR (htseq-counts)
countsTableRaw <- read.table("countsTOTALS_CodingGenes.tsv", header = TRUE)
caca  <- read.table("countsTOTALS.tsv", header = TRUE)

# Remove the total RNA.
countsTableNoto <- countsTableRaw[, 7:24]

#total analysis
countsTableTOT <- countsTableRaw[, 1:6]

# Generate the samples grouping.
groups <- as.factor(c("MonoH","MonoH","MonoH", "MonoL","MonoL","MonoL", "LightH","LightH","LightH", "LightL", "LightL","LightL", "HeavyH","HeavyH","HeavyH", "HeavyL","HeavyL","HeavyL"))
groupstot <- as.factor(c("TotH","TotH","TotH","TotL","TotL","TotL"))


# Filter lowly expressed genes according to our homemade function.
cpmTable <- lowExpression_filter(cpm(countsTableNoto), groups, tresh = 2.5)
cpmTableTot <- lowExpression_filter(cpm(countsTableTOT), groupstot, tresh = 2.5)
cpmall <- lowExpression_filter(cpm(countsTableRaw), groupsall, tresh = 2.5)


#rapresenting cpm single genes
groupsall <-  as.factor(c("Total","Total","Total","Total","Total","Total","Mono","Mono","Mono", "Mono","Mono","Mono", "Light","Light","Light", "Light", "Light","Light", "Heavy","Heavy","Heavy", "Heavy","Heavy","Heavy"));
treasAll <- as.factor(c(rep(c("High", "High", "High", "Low", "Low", "Low"), 4)))

barplot(colSums(cpmall), col = c("lightgreen", "lightgreen","lightgreen","lightblue", "lightblue", "lightblue"), main = "Library size CPM")

geneBoxplotCond(vm$E,"ENSG00000196262",groupsall,treasAll)

# Retrieve the original counts from the previous table.
countsTable <- countsTableNoto[rownames(cpmTable),]
countsTableTOT <- countsTableTOT[rownames(cpmTableTot),]

# Build the design matrix
designMat <- model.matrix(~0+groups);
colnames(designMat) <- levels(groups);
designMatTot <- model.matrix(~0+groupstot);

# Build the contrast matrix.
contr.matrix <- makeContrasts(
  MonoHvsMonoL = MonoH - MonoL,
  LightHvsLightL = LightH - LightL,
  HeavyHvsHeavyL = HeavyH - HeavyL,
  levels = colnames(designMat))

contr.matrixTot <- makeContrasts(
  TotHvsTotL = groupstotTotH - groupstotTotL,
  levels = colnames(designMatTot))

# Start the voom-limma analysis.
vm <- voom(countsTable, designMat)
fit <- lmFit(vm, designMat)
vfit <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(vfit, robust = TRUE) # Playing with the parameters of ebayes makes no difference.
ef <- decideTests(efit, lfc= 0.5)
summary(ef)
tfit <- treat(vfit,robust = TRUE, lfc =log2(1.1))
tf <- decideTests(tfit)
summary(tf)

#voom-limma analysis tot
vmtot <- voom(countsTableTOT, designMatTot)
fitTot <- lmFit(vmtot, designMatTot)
vfitTot <- contrasts.fit(fitTot, contrasts = contr.matrixTot)
tfittot <- treat(vfitTot,robust = TRUE, lfc = log2(1.15))
efitot <- eBayes(vfitTot, robust = TRUE)
tftot <- decideTests(tfittot)
eftot <- decideTests(efitot, lfc = 0.5)
summary(tftot)
summary(eftot)

plotMD(efitot, column = 1, status = eftot[,1], main = colnames(efitot)[1])
voomPCAtot <- PCA(t(vmtot$E))
barplot(colSums(countsTableTOT), col = c("lightgreen", "lightgreen","lightgreen","lightblue", "lightblue", "lightblue"), main = "Library size Raw counts")
barplot(colSums(vmtot$E), col = c("lightgreen", "lightgreen","lightgreen","lightblue", "lightblue", "lightblue"), main = "Normalised counts (CPM voom method)",ylab = "CPM")


# Venn diagram TFIT
vennDiagram(ef[,1:3], circle.col = c("turquoise", "salmon", "lightgreen"), main = "Venn with tFIT")


# MD plots TFIT
plotMD(efit, column = 1, status = ef[,1], main = colnames(efit)[1])
plotMD(efit, column = 2, status = ef[,2], main = colnames(efit)[2])
plotMD(efit, column = 3, status = ef[,3], main = colnames(efit)[3])


# Venn diagram EFIT
vennDiagram(ef[,1:3], circle.col = c("turquoise", "salmon", "lightgreen"), main = "Venn with EFIT")
vennDiagram(tf[,1:3], circle.col = c("turquoise", "salmon", "lightgreen"), main = "Venn with EFIT")

par(mfrow = c(1, 3))
# MD plots eFIT
plotMD(efit, column = 1, status = ef[,1], main = colnames(efit)[1])
plotMD(efit, column = 2, status = ef[,2], main = colnames(efit)[2])
plotMD(efit, column = 3, status = ef[,3], main = colnames(efit)[3])

# Series of plots for QC and visualisation of crap.
# Barplots of raw counts.
barplot(colSums(countsTable), col = c("lightgreen", "lightgreen","lightgreen","lightblue", "lightblue", "lightblue"), main = "Library size Raw counts")
barplot(colSums(vm$E), col = c("lightgreen", "lightgreen","lightgreen","lightblue", "lightblue", "lightblue"), main = "Normalised counts (CPM voom method)",ylab = "CPM")

# PCA on normalised counts.
voomPCA <- PCA(t(vm$E))
#voomPCA <- PCA(t(vm$E), graph = FALSE)
#fviz_pca_ind(voomPCA, habillage = ...., label="none", addEllipses = TRUE) # This produces a very nice plot so we keep it for publication only.
voomPCA <- PCA(t(vm$E))
# PCAs on each experiment.
PCA(t(vm$E[,1:6]))
PCA(t(vm$E[,7:12]))
PCA(t(vm$E[,13:18]))

# DEG tables.
DEheavy <- topTreat(efit, coef=1, n=Inf, p.value=0.05, lfc=1)
DElight <- topTreat(efit, coef=2, n=Inf, p.value=0.05, lfc=1)
DEmono <- topTreat(efit, coef=3, n=Inf, p.value=0.05, lfc=1)

# Retrive names of DEGs and pool them.
DEgenesefit <- union(rownames(DEheavy), rownames(DElight))
DEgenesefit <- union(DEgenestfit, rownames(DEmono))

selectef <- ef[DEgenesefit,]
vennDiagram(selectef[,1:2], circle.col = c("turquoise", "salmon", "lightgreen"), main = "Venn with EFIT")

# Subtract the DE genes from the cpm table.
genesDEcpm <- cpmTable[DEgenestfit,]

#Calculate means of groups
monoH <- rowMeans(genesDEcpm[,1:3])
monoL <- rowMeans(genesDEcpm[,4:6])
lightH <- rowMeans(genesDEcpm[,7:9])
lightL <- rowMeans(genesDEcpm[,10:12])
heavyH <- rowMeans(genesDEcpm[,13:15])
heavyL <- rowMeans(genesDEcpm[,16:18])

genesDEavg <- data.frame("monoH" = monoH, "monoL" = monoL, "lightH" = lightH,"lightL" = lightL,"heavyH"= heavyH, "heavyL" = heavyL, row.names = names(heavyH))

#Calculate logRatios over means
ratioH <- log2(genesDEavg$heavyH/genesDEavg$heavyL)
ratioL <- log2(genesDEavg$lightH/genesDEavg$lightL)
ratioM <- log2(genesDEavg$monoH/genesDEavg$monoL)

logRatiosDEGtable <- data.frame("Mono" = ratioM, "Light" = ratioL, "Heavy" = ratioH,  row.names = rownames(genesDEavg))

#Heatmap Hclustering
my_palette <- brewer.pal(n = 11, name = "RdYlGn")

heatmap.2(as.matrix(logRatiosDEGtable), col=my_palette, dendrogram = "row", cexCol = 1.5, cexRow = 0.0001, key.title = NA, keysize = 0.9, main = "DEGs logRatio H/l Clustering", ylab = "Genes")

#determine number of cluster by plotting sum of squares whitnin groups
wss <- (nrow(logRatiosDEGtable)-1)*sum(apply(logRatiosDEGtable,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(logRatiosDEGtable,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


# Perform k-means clustering k-means
k <- kmeans(as.matrix(logRatiosDEGtable), 9)


# Append id and cluster
dfc <- cbind(logRatiosDEGtable, id=seq(nrow(logRatiosDEGtable)), cluster=k$cluster)

# Add idsort, the id number ordered by cluster
dfc$idsort <- dfc$id[order(dfc$cluster)]
dfc$idsort <- order(dfc$idsort)

clusterCols <- as.character(sort(k$cluster))
#clusterCols <- str_replace_all(clusterCols, "1", "grey") #to change colors

heatmap(as.matrix(logRatiosDEGtable)[order(k$cluster),], Rowv=NA,Colv=NA, scale="none",labRow=NA, cexCol = 1.5, col=my_palette, RowSideColors = clusterCols, ylab= "Genes LogRatio", main= "k-means(7) clustering of H/L LogRatios")


cluster1 <- rownames(dfc[which(dfc$cluster == 1),])
cluster2 <- rownames(dfc[which(dfc$cluster == 2),])
cluster3 <- rownames(dfc[which(dfc$cluster == 3),])
cluster4 <- rownames(dfc[which(dfc$cluster == 4),])
cluster5 <- rownames(dfc[which(dfc$cluster == 5),])
cluster6 <- rownames(dfc[which(dfc$cluster == 6),])


write(cluster1, "cluster1_Kmeans.txt")
write(cluster2, "cluster2_Kmeans.txt")
write(cluster3, "cluster3_Kmeans.txt")
write(cluster4, "cluster4_Kmeans.txt")
write(cluster5, "cluster5_Kmeans.txt")
write(cluster6, "cluster6_Kmeans.txt")
