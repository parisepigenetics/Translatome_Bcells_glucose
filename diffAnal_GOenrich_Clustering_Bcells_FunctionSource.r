# Source file with usefull finctions for the diffAnal_EnrichAnal_ClustAnal of the Beta-cells polysome profile.
# cbouyio, May 2019, UMR 7216

# Functions -------------------------------------------------------------------

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
  # Function to filter lowly expressed genes.
  # e      : raw counts data.frame (or cpm, or tpm matrix).
  # groups : factor designating the grouping(s) (conditions, treatments etc.) it MUSt be of equal length to the columns of e and its levels must represent the different groups (conditions).
  # thres  : the threshold of the *mean* between groups.
  # samples: the minimum number of groups (i.e. conditions) that we want the mean to be higher than the threshold 'thres'.
  # coefVar: The coefficient of variation threshold, within replicates, used to remove "noisy" measured genes.
  filteredDF <- data.frame()
  rows <- vector()
  for (i in 1:nrow(e)) {
    row <- as.numeric(e[i,])
    means <- aggregate(row~groups, FUN = mean)$row
    sds <- aggregate(row~groups, FUN = sd)$row
    cvs <- sds/means
    lg <- length(levels(groups))
    # This condition filters for "noisy" genes and for lowly expressed genes.
    if ( (sum(cvs <= coefVar) == lg) & (sum(means >= thres) >= samples) ) { # The first part of this condition checks for the coefficient of variability in ALL groups by grouping sums in the number of samples and the second part is checking for the experiment threshold we specify.
      filteredDF <- rbind(filteredDF, row)
      rows <- append(rows, rownames(e)[i])
    }
  }
  colnames(filteredDF) <- colnames(e)
  rownames(filteredDF) <- rows
  return(filteredDF)
}


plot_semiSupervised_clust <- function(data, k, method, scale = FALSE, title = "", ...){
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
  # Title
  if (title == "") {
    ti = paste(method, " clustering of, ", deparse(substitute(data)), " scale: ", as.character(scale))
  } else {
    ti = title
  }
  # Plotting
  heatmap(as.matrix(data)[order(clustRes$cluster),], Rowv = NA, Colv = NA, scale = "none", labRow = NA, cexCol = 1.5, col = my_palette, RowSideColors = clusterCols, ylab = "Genes", main = ti)
  invisible(list(res = clustRes, df = dfcall))
}


plot_unSupervised_clust <- function(data, method, scale = FALSE, title = "", ...){
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
  # Title
  if (title == "") {
    ti = paste(method, " clustering of, ", deparse(substitute(data)), " scale: ", as.character(scale))
  } else {
    ti = title
  }
  # Plotting
  heatmap(as.matrix(data)[order(clustRes$classification),], Rowv = NA, Colv = NA, scale = "none", labRow = NA, cexCol = 1.5, col = my_palette, RowSideColors = clusterCols, ylab = "Genes", main = ti)
  invisible(list(res = clustRes, df = dfcall))
}


# Low lwvwel function to put the counts on a boxplot.
n_fun <- function(x){
  return(data.frame(y = max(x), label = length(x)))
}
