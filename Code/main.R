library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)

### Determining the Data Series and Platform
series = "GSE48558"
platform = "GPL6244"

### Downloading the Dataset from GEO
gSet = GEOquery::getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gSet) > 1) idx = grep("GPL6244", attr(gSet, "names")) else idx <- 1
gSet = gSet[[idx]]
  
### Printing the Head data
head(gSet)

### Grouping Samples based on their Source Name
gr = c('AMLPatient', 'AMLPatient', 'AMLPatient', 'AMLPatient', 'AMLPatient', 'AMLPatient', 'AMLPatient', 'AMLPatient', 'AMLPatient', 'AMLPatient', 'AMLPatient', 'AMLPatient', 'AMLPatient', 'BALLCellLine', 'BALLCellLine', 'BALLCellLine', 'BALLCellLine', 'TALLCellLine', 'TALLCellLine', 'BALLCellLine', 'TALLCellLine', 'BALLCellLine', 'BALLCellLine', 'TALLCellLine', 'BALLCellLine', 'BALLCellLine', 'TALLCellLine', 'TALLCellLine', 'BALLCellLine', 'BALLCellLine', 'BALLCellLine', 'BALLCellLine', 'BALLCellLine', 'TALLCellLine', 'TALLCellLine', 'BALLCellLine', 'TALLCellLine', 'BALLPatient', 'TALLCellLine', 'AMLCellLine', 'Granulocytes', 'BALLPatient', 'TALLCellLine', 'AMLCellLine', 'Granulocytes', 'BALLPatient', 'TALLCellLine', 'AMLCellLine', 'BALLPatient', 'AMLCellLine', 'AMLCellLine', 'BALLPatient', 'AMLCellLine', 'AMLCellLine', 'BALLPatient', 'AMLCellLine', 'AMLCellLine', 'BALLPatient', 'BALLCellLine', 'AMLCellLine', 'BALLPatient', 'BALLCellLine', 'AMLCellLine', 'BALLPatient', 'BALLCellLine', 'AMLCellLine', 'BALLPatient', 'BALLCellLine', 'BCells', 'BALLCellLine', 'TCells', 'AMLCellLine', 'BALLPatient', 'BALLCellLine', 'Granulocytes', 'BALLCellLine', 'Granulocytes', 'Monocytes', 'Monocytes', 'BCells', 'BALLCellLine', 'TCells', 'AMLCellLine', 'BALLCellLine', 'TCells', 'TCells', 'AMLCellLine', 'BALLCellLine', 'TCells', 'TCells', 'AMLCellLine', 'BCells', 'BALLCellLine', 'TCells', 'AMLCellLine', 'BCells', 'BALLCellLine', 'TCells', 'AMLCellLine', 'CD34', 'TALLCellLine', 'TALLPatient', 'AMLCellLine', 'CD34', 'TALLCellLine', 'TALLPatient', 'AMLCellLine', 'CD34', 'TALLCellLine', 'TALLPatient', 'AMLCellLine', 'TALLPatient', 'BALLPatient', 'BALLPatient', 'BALLPatient', 'BALLPatient', 'BALLPatient', 'BALLPatient', 'BALLPatient', 'BALLPatient', 'BALLPatient', 'TALLPatient', 'BALLPatient', 'BALLPatient', 'BALLPatient', 'BALLPatient', 'BALLPatient', 'BALLPatient', 'BALLPatient', 'TALLPatient', 'TALLPatient', 'TALLPatient', 'TALLPatient', 'TALLPatient', 'TALLPatient', 'TALLPatient', 'TALLPatient', 'Granulocytes', 'Granulocytes', 'Granulocytes', 'Granulocytes', 'Granulocytes', 'Granulocytes', 'Granulocytes', 'AMLPatient', 'AMLPatient', 'TCells', 'AMLPatient', 'AMLPatient', 'AMLPatient', 'BCells', 'BCells', 'BCells', 'BCells', 'BCells', 'BCells', 'BCells', 'TCells', 'Monocytes', 'Monocytes', 'Monocytes', 'Monocytes', 'Granulocytes', 'TCells', 'TCells', 'TCells', 'TCells', 'TCells', 'TCells', 'TCells')
print(gr)

### Expression Matrix
exMat <- Biobase::exprs(gSet)

print(gSet)
print(dim(exMat))

### Plotting the Expression Matrix
pdf("D:/Main/Uni/TERM5/Bio/Project/Result/boxplot.pdf", width = 15, height = 15)
boxplot(exMat)
dev.off()

### Plotting the Heat Map, where differences between Genes are indicated by Color
pdf("D:/Main/Uni/TERM5/Bio/Project/Result/CorHeatmapGr.pdf", width = 20, height = 20)
pheatmap(cor(exMat), labels_row = gr, labels_col = gr, border_color = NA)
dev.off()

### Implementing and Plotting Principle Component Analysis
pc = prcomp(exMat)
pdf("D:/Main/Uni/TERM5/Bio/Project/Result/First Principal Component.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

### Converting Principle Component Analysis to Scaled one
exMat.scale = t(scale(t(exMat), scale = FALSE))
mean(exMat.scale[1,])
scaled_pc = prcomp(exMat.scale)
pdf("D:/Main/Uni/TERM5/Bio/Project/Scaled PC.pdf")
plot(scaled_pc)
plot(scaled_pc$x[,1:2])
dev.off()

### Dimension Reduction using PCA and Plotting on a 2-D plane
pcr = data.frame(pc$rotation[,1:3], Group = gr)
pdf("D:/Main/Uni/TERM5/Bio/Project/Result/PCA Samples.pdf")
ggplot(pcr, aes(x = PC1, y = PC2, color = Group)) + geom_point(size = 3) + theme_bw()
dev.off()


### Differential Expression Analysis
gsms <- paste0("0000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXX1XXX1XXXXX",
               "XXXXXXXXXXXXXXXXXX1X1XXX1X1111X1XX11XX11X1X1X1X1X1",
               "XXX1XXX1XXXXXXXXXXXXXXXXXXXXXXXXXXXXX1111111001000",
               "11111111111111111111")
sml <- strsplit(gsms, split="")[[1]]
sel <- which(sml != "X")
sml <- sml[sel]
gSet <- gSet[ ,sel]

gr <- factor(sml)
groups <- make.names(c("AML","Healthy"))
levels(gr) <- groups
gSet$group <- gr
design <- model.matrix(~group + 0, gSet)
colnames(design) <- levels(gr)
fit <- lmFit(gSet, design)  # fit linear model


cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("Gene.symbol","Gene.ID","adj.P.Val","logFC"))
write.table(tT, file="D:/Main/Uni/TERM5/Bio/Project/Result/AML_Normal.txt", row.names=F, sep="\t", quote=F)


aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(as.character((strsplit2(aml.up$Gene.symbol, "///"))))
write.table(aml.up.genes, file="D:/Main/Uni/TERM5/Bio/Project/Result/AML_Normal_Up.txt", row.names=F,
            col.names = F,sep="\t", quote=F)
aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character((strsplit2(aml.down$Gene.symbol, "///"))))
write.table(aml.down.genes, file="D:/Main/Uni/TERM5/Bio/Project/Result/AML_Normal_Down.txt", row.names=F,
            col.names = F,sep="\t", quote=F)
