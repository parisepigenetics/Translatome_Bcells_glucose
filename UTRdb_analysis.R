library(tidyverse)
library(ggplot2)
library(gplots)
library(pheatmap)

#read table, befor remove coordinates in feature column[]
utr_5 <- read.table("utrScan_5utr_results.txt", sep = ":", header = TRUE)
utr_3 <- read.table("utrScan_3utr_results.txt", sep = ":", header = TRUE)


utr_5 <- as.data.frame(utr_5[,1:2])
utr_3 <- as.data.frame(utr_3[,1:2])

utr5_1 <- utr_5 %>% count(trascript_ID, feature, sort = TRUE, .drop= FALSE)
utr3_1 <- utr_3 %>% count(trascript_ID, feature, sort = TRUE, .drop= FALSE)

utr5_1 <- as.data.frame(utr5_1)
utr3_1 <- as.data.frame(utr3_1)

utr5_matrix_features <- utr5_1 %>%  spread(key = feature, value = n)
utr3_matrix_features <- utr3_1 %>%  spread(key = feature, value = n)

rownames(utr5_matrix_features) <- utr5_matrix_features$trascript_ID
utr5_matrix_features$trascript_ID <- NULL

rownames(utr3_matrix_features) <- utr3_matrix_features$trascript_ID
utr3_matrix_features$trascript_ID <- NULL
utr3_matrix_features$uORF <- NULL

library(RColorBrewer)
coul = colorRampPalette(brewer.pal(8, "Blues"))(13)

heatmap(as.matrix(utr5_matrix_features), col = coul)
heatmap(as.matrix(utr3_matrix_features))
pheatmap(as.matrix(utr5_matrix_features))


