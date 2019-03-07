library(ggplot2)
library(DESeq2)
library(reshape)
library(VennDiagram)
library(pheatmap)
library(limma)
library(Glimma)
library(edgeR)
library(FactoMineR)
library(factoextra)
library(gplots)
library(RColorBrewer)
library(ggdendro)

###Functions

# Function to filter form lowly expressed genes.
# Function to filter lowly expressed genes.
lowExpression_filter <- function(e, groups, thres = 3, coefVar = 0.5, samples = 2, ...){
  # e : raw counts data.frame (or cpm, or better tpm matrix).
  # groups : factor designating the grouping(s) (conditions, treatments etc.)
  #it MUSt be of equal lenth to the columns of e and its levels must be the  
  #different groups.
  # thres : the threshold of the *mean* between groups.
  # samples : the minimum number of sample groups that we want threshold
  #thres to be higher.
  filteredDF <- data.frame()
  rows <- vector()
  for (i in 1:nrow(e)) {
    row <- as.numeric(e[i,])
    means <- aggregate(row~groups, FUN = mean)$row
    sds <- aggregate(row~groups, FUN = sd)$row
    cvs = sds/means
    ll <- length(levels(groups))
    if ( (sum(cvs <= coefVar) == ll) & (sum(means >= thres) >=samples) ) { # This condition checks for the coefficient of variability in ALL groups and the grouping sums in the number of samples.
      filteredDF <- rbind(filteredDF, row)
      rows <- append(rows, rownames(e)[i])
    }
  }
  colnames(filteredDF) <- colnames(e)
  rownames(filteredDF) <- rows
  return(filteredDF)
}

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

#START
DEGenes <- scan("final_DEgenes.txt", what = "character")

tmpmeans <-  read.table("polysomeProfile_TPM_proteinCoding_CleanedReliable.csv",sep = ";", header = T)

tpmall1 <- read.table("polysomeProfile_TPM_proteinCoding.csv", header = TRUE)
tpmall2 <- tpmall1[,1:18]
tpmall3 <- tpmall1[,1:6]

countsTableRaw <- read.table("countsTOTALS_CodingGenes.tsv", header = TRUE)
countsTableRaw1 <- countsTableRaw[,7:24]


countsTableRaw <- countsTableRaw[c("MonoH1", "MonoH2", "MonoH3", "MonoL1", "MonoL2", "MonoL3",
                                   "LightH1", "LightH2", "LightH3", "LightL1", "LightL2", "LightL3",
                                   "HeavyH1", "HeavyH2", "HeavyH3", "HeavyL1", "HeavyL2", "HeavyL3",
                                   "TotH1", "TotH2", "TotH3", "TotL1", "TotL2", "TotL3")];

groupsall <- factor(c("Mono","Mono","Mono","Mono","Mono","Mono",
                      "Light","Light","Light","Light","Light","Light",
                      "Heavy","Heavy","Heavy","Heavy","Heavy","Heavy",
                      "Total","Total","Total","Total","Total","Total"),
                    levels = c("Mono", "Light", "Heavy", "Total"));

treasAll <- as.factor(c(rep(c("High", "High", "High", "Low", "Low", "Low"), 4)));

# Generate the samples grouping.
groupsformatrix <- factor(c("MonoH","MonoH","MonoH","MonoL","MonoL","MonoL","LightH","LightH",
                          "LightH","LightL","LightL","LightL","HeavyH","HeavyH","HeavyH","HeavyL","HeavyL",
                          "HeavyL","TotalH","TotalH","TotalH","TotalL","TotalL","TotalL"),  
                          levels = c("TotalH","TotalL","MonoH","MonoL", "LightH","LightL", "HeavyH","HeavyL"))

groupsformatrix_notot <- factor(c("MonoH","MonoH","MonoH","MonoL","MonoL","MonoL","LightH","LightH",
                            "LightH","LightL","LightL","LightL","HeavyH","HeavyH","HeavyH","HeavyL","HeavyL",
                            "HeavyL"),  
                          levels = c("MonoH","MonoL", "LightH","LightL", "HeavyH","HeavyL"))

# Filter lowly expressed genes according to our homemade function.
cpmtmp1 <- cpm(countsTableRaw)
cpmtmp <- cpmtmp1[apply(cpmtmp1>1, 1, all),]
cpmall <- lowExpression_filter(cpmtmp, groupsformatrix, thres = 2.5, samples = 6 )
cpmall <- lowExpression_filter(cpmtmp, groupsformatrix_notot, thres = 4, samples = 6)

dim(cpmall1)
dim(cpmall)
summary(cpmall)

countsTableall <- countsTableRaw[rownames(cpmall),]
countsTable_noTOT <- countsTableall[,7:24]
countsTable_noTOT <- countsTableRaw1[rownames(cpmall),]

tpmtrim <- tpmall2[rownames(cpmall),]
tpmtrimT <- tpmall3[rownames(cpmall),]

#design
designMatall <- model.matrix(~0+groupsformatrix);
colnames(designMatall) <- levels(groupsformatrix);

designMat_notot <- model.matrix(~0+groupsformatrix_notot);
colnames(designMat_notot) <- levels(groupsformatrix_notot);

# Build the contrast matrix.
contr.matrixAll <- makeContrasts(
  TotHvsTotL = TotalH - TotalL,
  monoHvsmonoL = MonoH - MonoL,
  lightHvslightL = LightH - LightL,
  heavyLvsheavyH = HeavyH - HeavyL,
  levels = colnames(designMatall))

contr.matrix_notot <- makeContrasts(
   monoHvsmonoL = MonoH - MonoL,
  lightHvslightL = LightH - LightL,
  heavyLvsheavyH = HeavyH - HeavyL,
  levels = colnames(designMat_notot))

# Start the voom-limma analysis.
###NO total
vmnotot <- voom(countsTable_noTOT, designMat_notot, plot = TRUE)
fitnotot <- lmFit(vmnotot, designMat_notot)
vfitnotot <- contrasts.fit(fitnotot, contrasts = contr.matrix_notot)
efitnotot <- eBayes(vfitnotot, robust = TRUE) # Playing with the parameters of ebayes makes no difference.
efnotot <- decideTests(efitnotot, lfc = 0.5, p.value = 0.05)
summary(efnotot)
plotSA(efitnotot)
tfitnotot <- treat(vfitnotot,robust = TRUE, lfc= log2(1.15))
ttnotot <- decideTests(tfitnotot)
summary(ttnotot)



vmall <- voom(countsTableall, designMatall, plot = TRUE)
fitall <- lmFit(vmall, designMatall)
vfitall <- contrasts.fit(fitall, contrasts = contr.matrixAll)
efitall <- eBayes(vfitall) # Playing with the parameters of ebayes makes no difference.
efall <- decideTests(efitall, lfc = 0.5, p.value = 0.05)
summary(efall)
plotSA(efitall)
tfitall <- treat(vfitall,robust = TRUE, lfc= log2(1.15))
ttall <- decideTests(tfitall)
summary(ttall)

a<-topTreat(tfitall, coef=4, n=Inf, p.value=0.05, lfc= 0.5)


plotMD(efitall, column = 1, status = efall[,1], main = colnames(efitall)[1])
plotMD(efitall, column = 2, status = efall[,2], main = colnames(efitall)[2])
plotMD(efitall, column = 3, status = efall[,3], main = colnames(efitall)[3])
plotMD(efitall, column = 4, status = efall[,4], main = colnames(efitall)[4])

vennDiagram(efall[,2:4], circle.col = c("turquoise", "salmon", "lightgreen"), main = "Venn with TFIT")

vennDiagram(efnotot, circle.col = c("turquoise", "salmon", "lightgreen"), main = "Venn with TFIT")

barplot(colSums(countsTableall), col = c("lightgreen", "lightgreen","lightgreen","lightblue", "lightblue", "lightblue"), main = "Library size Raw counts")
barplot(colSums(vmall$E), col = c("lightgreen", "lightgreen","lightgreen","lightblue", "lightblue", "lightblue"), main = "Normalised counts (CPM voom method)",ylab = "CPM")

voomPCA <- PCA(t(vmall$E))
voomPCA <- PCA(t(vmnotot$E))
plotMDS(vmnotot$E)

dist_mat <- get_dist(t(tpmtrim), method = 'euclidian')
hclust_avg <- hclust(dist_mat, method = "ward.D2")
fviz_dend(hclust_avg, k = 6, cex = 1.25, k_colors = c("dodgerblue4","firebrick","firebrick", "dodgerblue4","dodgerblue4","firebrick"),
                                                color_labels_by_k = TRUE, rect = TRUE, lower_rect = -6000,
                                                main = "Hierarchical Clustering Sequenced Samples")

dist_mat <- get_dist(t(tpmtrimT), method = 'euclidian')
hclust_avg <- hclust(dist_mat, method = "ward.D2")
fviz_dend(hclust_avg, k = 2, cex = 1.25, k_colors = c("dodgerblue4","firebrick"),
          color_labels_by_k = TRUE, rect = TRUE, lower_rect = -1000,
          main = "Hierarchical Clustering Transcriptome Samples")


PCA(t(vmall$E[,1:6]))
PCA(t(vmall$E[,7:12]))
PCA(t(vmall$E[,13:18]))
PCA(t(vmall$E[,19:24]))

DEheavyall <- topTable(efitall, coef=4, n=Inf, p.value=0.05, lfc = 0.3)
DElightall <- topTable(efitall, coef=3, n=Inf, p.value=0.05, lfc = 0.3)
DEmonoall <- topTable(efitall, coef=2, n=Inf, p.value=0.05, lfc = 0.3)
DEtotal <- topTable(efitall, coef=1, n=Inf, p.value=0.05, lfc = 0.3)

###
logFCMono <- data.frame ("Monosomes" =DEmonoall$logFC, row.names = rownames(DEmonoall))
logFClight<- data.frame ("Light_poly" =DElightall$logFC, row.names = rownames(DElightall))
logFCheavy <- data.frame ("Heavy_poly" =DEheavyall$logFC, row.names = rownames(DEheavyall))
logFCtot <- data.frame ("Total" =DEtotal$logFC, row.names = rownames(DEtotal))


logFCall <- data.frame("Total"= logFCtot,"Monosomes"= logFCMono,"Light_poly"= logFClight, "Heavy_poly"= logFCheavy)

logFCall1 <- reshape2::melt(logFCall, id.vars = NULL)


logFCall1$LogFC <- 2^logFCall1$LogFC

p <- ggplot(logFCall1, aes(x = Condition, y = LogFC, fill= Condition)) + geom_violin() + geom_boxplot(width=0.1)


####

DEgenestfitall <- union(rownames(DEheavyall), rownames(DElightall))
DEgenestfitall <- union(DEgenestfitall, rownames(DEmonoall))

selectetall <- ttall[DEgenestfitall,]
vennDiagram(selectetall[,2:4], circle.col = c("turquoise", "salmon", "lightgreen"), main = "Venn with t-FIT")


genesDEcpmall <- cpmall[DEgenestfitall,]
genesDEtpmall <- tpmall1[DEGenes,]
logDEmono <- logFCMono[row.names(logFCMono),]

#Calculate means of groups
monoHall <- rowMeans(genesDEcpmall[,1:3])
monoLall <- rowMeans(genesDEcpmall[,4:6])
lightHall <- rowMeans(genesDEcpmall[,7:9])
lightLall <- rowMeans(genesDEcpmall[,10:12])
heavyHall <- rowMeans(genesDEcpmall[,13:15])
heavyLall <- rowMeans(genesDEcpmall[,16:18])

monoHall1 <- rowMeans(genesDEtpmall[,1:3])
monoLall1 <- rowMeans(genesDEtpmall[,4:6])
lightHall1 <- rowMeans(genesDEtpmall[,7:9])
lightLall1 <- rowMeans(genesDEtpmall[,10:12])
heavyHall1 <- rowMeans(genesDEtpmall[,13:15])
heavyLall1 <- rowMeans(genesDEtpmall[,16:18])

genesDEavgall <- data.frame("monoHall" = monoHall, "monoLall" = monoLall, "lightHall" = lightHall,"lightLall" = lightLall,"heavyHall"= heavyHall, "heavyLall" = heavyLall, row.names = names(heavyHall))

genesDEavgall1 <- data.frame("monoHall1" = monoHall1, "monoLall1" = monoLall1, "lightHall1" = lightHall1,"lightLall1" = lightLall1,"heavyHall1"= heavyHall1, "heavyLall" = heavyLall1, row.names = names(heavyHall1))


#Calculate logRatios over means
ratioHall <- log2(genesDEavgall$heavyHall/genesDEavgall$heavyLall)
ratioLall <- log2(genesDEavgall$lightHall/genesDEavgall$lightLall)
ratioMall <- log2(genesDEavgall$monoHall/genesDEavgall$monoLall)

ratioHall1 <- log2(genesDEavgall1$heavyHall/genesDEavgall1$heavyLall)
ratioLall1 <- log2(genesDEavgall1$lightHall/genesDEavgall1$lightLall)
ratioMall1 <- log2(genesDEavgall1$monoHall/genesDEavgall1$monoLall)

logRatiosDEGtableall <- data.frame("Mono" = ratioMall, "Light" = ratioLall, "Heavy" = ratioHall,  row.names = rownames(genesDEavgall))

logRatiosDEGtableall1 <- data.frame("Mono" = ratioMall1, "Light" = ratioLall1, "Heavy" = ratioHall1,  row.names = rownames(genesDEavgall1))

logRatiosDEGtableall2 <- logRatiosDEGtableall1[,2:3]

#Heatmap Hclustering
my_palette <- brewer.pal(n = 10, name = "RdYlGn")

heatmap.2(as.matrix(logRatiosDEGtableall), col=my_palette, dendrogram = "row", cexCol = 1.5, cexRow = 0.0001, key.title = NA, keysize = 0.9, main = "DEGs logRatio H/l Clustering", ylab = "Genes")


heatmap.2(as.matrix(logRatiosDEGtableall1) ,Colv = "as-is" , col=my_palette, dendrogram = "row", cexCol = 1.5, cexRow = 0.0001, key.title = NA, keysize = 0.9, main = "DEGs logRatio H/l Clustering", ylab = "Genes")

heatmap.2(as.matrix(logRatiosDEGtableall2) ,Colv = "as-is" , col=my_palette, dendrogram = "row", cexCol = 1.5, cexRow = 0.0001, key.title = NA, keysize = 0.9, main = "DEGs logRatio H/l Clustering", ylab = "Genes")


#determine number of cluster by plotting sum of squares whitnin groups
wss <- (nrow(logRatiosDEGtableall)-1)*sum(apply(logRatiosDEGtableall,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(logRatiosDEGtableall,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


wss <- (nrow(logRatiosDEGtableall)-1)*sum(apply(logRatiosDEGtableall,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(logRatiosDEGtableall1,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


# Perform k-means clustering k-means
k <- kmeans(as.matrix(logRatiosDEGtableall), 7)
k1 <- kmeans(as.matrix(logRatiosDEGtableall1), 3)

fviz_cluster(k1, logRatiosDEGtableall1, geom = "point")
fviz_cluster(k, logRatiosDEGtableall, geom = "point")

# Append id and cluster
dfcall <- cbind(logRatiosDEGtableall, id=seq(nrow(logRatiosDEGtableall)), cluster=k$cluster)

dfcall1 <- cbind(logRatiosDEGtableall1, id=seq(nrow(logRatiosDEGtableall1)), cluster=k1$cluster)

# Add idsort, the id number ordered by cluster
dfcall$idsort <- dfcall$id[order(dfcall$cluster)]
dfcall$idsort <- order(dfcall$idsort)

dfcall1$idsort <- dfcall1$id[order(dfcall1$cluster)]
dfcall1$idsort <- order(dfcall1$idsort)

clusterCols <- as.character(sort(k$cluster))

clusterCols1 <- as.character(sort(k1$cluster))

#clusterCols <- str_replace_all(clusterCols, "1", "grey") #to change colors

heatmap(as.matrix(logRatiosDEGtableall)[order(k$cluster),], Rowv=NA,Colv=NA, scale="none",labRow=NA, cexCol = 1.5,
        col=my_palette, RowSideColors = clusterCols, ylab= "Genes LogRatio",
        main= "k-means(7) clustering of H/L LogRatios")

heatmap(as.matrix(logRatiosDEGtableall1)[order(k1$cluster),], Rowv=NA, Colv=NA, scale="none",labRow=NA, cexCol = 1.5,
        col=my_palette, RowSideColors = clusterCols1, ylab= "Genes LogRatio",
        main= "k-means() clustering of H/L LogRatios")

pheatmap(as.matrix(logRatiosDEGtableall1)[order(k1$cluster),], Rowv=NA, Colv=NA, scale="none",labRow=NA, cexCol = 1.5,
         col=my_palette, RowSideColors = clusterCols1, ylab= "Genes LogRatio", show_rownames = FALSE, cluster_cols = F,
         main= "k-means() clustering of H/L LogRatios")

cluster1 <- rownames(dfcall[which(dfcall$cluster == 1),])
cluster2 <- rownames(dfcall[which(dfcall$cluster == 2),])
cluster3 <- rownames(dfcall[which(dfcall$cluster == 3),])
cluster4 <- rownames(dfcall[which(dfcall$cluster == 4),])
cluster5 <- rownames(dfcall[which(dfcall$cluster == 5),])
cluster6 <- rownames(dfcall[which(dfcall$cluster == 6),])
cluster7 <- rownames(dfcall[which(dfcall$cluster == 7),])
cluster8 <- rownames(dfcall[which(dfcall$cluster == 8),])



vmcpmall <- data.frame(vmall$E)


#insulin
geneBoxplotCond(vmcpmall,"ENSG00000254647",groupsall,treasAll)
#PPIA
geneBoxplotCond(vmcpmall,"ENSG00000196262",groupsall,treasAll)
#ATF4
geneBoxplotCond(vmcpmall,"ENSG00000128272",groupsall,treasAll)
#CHOP
geneBoxplotCond(vmcpmall,"ENSG00000175197",groupsall,treasAll)
#JUN ENSG00000177606
geneBoxplotCond(vmcpmall,"ENSG00000177606",groupsall,treasAll)
#FOS? ENSG00000170345
geneBoxplotCond(vmcpmall,"ENSG00000170345",groupsall,treasAll)
#EGR1 ENSG00000120738
geneBoxplotCond(vmcpmall,"ENSG00000120738",groupsall,treasAll)
#RACK1
geneBoxplotCond(vmcpmall,"ENSG00000204628",groupsall,treasAll)
#nucleophospmin
geneBoxplotCond(cpmall,"ENSG00000181163",groupsall,treasAll)
#PABP
geneBoxplotCond(cpmall,"ENSG00000070756",groupsall,treasAll)
#HNRNPA1



cluster3b <- vmcpmall[cluster3,]



for (i in rownames(cac) ) {g <- geneBoxplotCond(tpmall1,i, groupsall,treasAll); show(g)}
show(g)
