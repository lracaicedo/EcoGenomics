# Laura Caicedo-Quiroga
# Script for HW#2

library("DESeq2")

library("ggplot2")

countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)

#subset intertidal and subtidal
colDataINT <-subset(colData, colData$location=="int")
colDataSUB <-subset(colData, colData$location=="sub")

countDataINT <-countData[, which(colnames(countData) %in% row.names(colDataINT))]
countDataSUB <-countData[, -which(colnames(countData) %in% row.names(colDataINT))]
dim(countDataINT)
dim(countDataSUB)

#################### MODEL NUMBER 1: TEST EFFECT OF HEALTH 
#First for Intertidal

ddsINT <- DESeqDataSetFromMatrix(countData = countDataINT, colData = colDataINT, design = ~ health)

dim(ddsINT)
ddsINT <- ddsINT[ rowSums(counts(ddsINT)) > 100, ]
dim(ddsINT)

colData(ddsINT)$health <- factor(colData(ddsINT)$health, levels=c("H","S")) 
#sets that "healthy is the reference

ddsINT <- DESeq(ddsINT) 

resINT <- results(ddsINT)
resINT <- resINT[order(resINT$padj),]
head(resINT)
summary(resINT)

#Intertidal
plotMA(resINT, main="DESeq2INT", ylim=c(-2,2))

## Check out one of the genes to see if it's behaving as expected....
d <- plotCounts(ddsINT, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","day")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + ylim(0,500)
p

## Check out one of the genes to see interaction between score, health and expression....
d <- plotCounts(ddsINT, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","score")), returnData=TRUE)
d
p <- ggplot(d, aes(x= score, y=count, shape = health)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p
p <- ggplot(d, aes(x=score, y=count, color=health, group=health)) 
p <- p +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10() + ggtitle("Intertidal")
p

#Model 2 for Subtidal

ddsSUB <- DESeqDataSetFromMatrix(countData = countDataSUB, colData = colDataSUB, design = ~ health)

dim(ddsSUB)
ddsSUB <- ddsSUB[ rowSums(counts(ddsSUB)) > 100, ]
dim(ddsSUB)

colData(ddsSUB)$health <- factor(colData(ddsSUB)$health, levels=c("H","S")) 
#sets that "healthy is the reference

ddsSUB <- DESeq(ddsSUB) 

resSUB <- results(ddsSUB)
resSUB <- resSUB[order(resSUB$padj),]
head(resSUB)
summary(resSUB)

#Subtidal
plotMA(resSUB, main="DESeq2SUB", ylim=c(-2,2))

## Check out one of the genes to see if it's behaving as expected....
d <- plotCounts(ddsSUB, gene="TRINITY_DN42073_c0_g1_TRINITY_DN42073_c0_g1_i1_g.12173_m.12173", intgroup=(c("health","day")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + ylim(0,500)
p

## Check out one of the genes to see interaction between score, health and expression....
d <- plotCounts(ddsSUB, gene="TRINITY_DN42073_c0_g1_TRINITY_DN42073_c0_g1_i1_g.12173_m.12173", intgroup=(c("health","score")), returnData=TRUE)
d
p <- ggplot(d, aes(x= score, y=count, shape = health)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p
p <- ggplot(d, aes(x=score, y=count, color=health, group=health)) 
p <- p +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10() + ggtitle("Subtidal")
p

### Model 3 accounting for location
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ location + health)
dim(dds)
dds <- dds[ rowSums(counts(dds)) > 100, ]
dim(dds)

colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S")) # (factor)sets that "healthy is the reference

dds <- DESeq(dds) 

res <- results(dds)
res <- res[order(res$padj),]
head(res)
summary(res)

# ALl accounting for location
# MA-Plot
plotMA(res, main="DESeq2ALL", ylim=c(-3,3.5))


## Check out one of the genes to see if it's behaving as expected....
d <- plotCounts(dds, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","day","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count, shape = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) + ylim(0,500)
p
## Check out one of the genes to see interaction between score, health and expression....
d <- plotCounts(dds, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","score","location")), returnData=TRUE)
d
p <- ggplot(d, aes(x= score, y=count, shape = health, color = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p
p <- ggplot(d, aes(x=score, y=count, color=health, group=health)) 
p <- p +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10() + ggtitle("All")
p



############## PCA plots
# For all
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("score"))
plotPCA(vsd, intgroup=c("health"))
plotPCA(vsd, intgroup=c("day"))
plotPCA(vsd, intgroup=c("location"))
plotPCA(vsd, intgroup=c("health","location"))


# For Intertidal
vsdINT <- varianceStabilizingTransformation(ddsINT, blind=FALSE)

plotPCA(vsdINT, intgroup=c("score"))
pca2<-plotPCA(vsdINT, intgroup=c("health"))
pca2 + ggtitle("Intertidal")
plotPCA(vsdINT, intgroup=c("day"))



# For Subtidal
vsdSUB <- varianceStabilizingTransformation(ddsSUB, blind=FALSE)

plotPCA(vsdSUB, intgroup=c("score"))
plotPCA(vsdSUB, intgroup=c("health"))
plotPCA(vsdSUB, intgroup=c("day"))



