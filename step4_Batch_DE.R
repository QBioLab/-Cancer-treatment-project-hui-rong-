library("sva")
library("DESeq2")
library("RColorBrewer")
library("IHW")
library("clusterProfiler")
library("org.Mm.eg.db")
library("enrichplot")
library("stringr")
library("pheatmap")

is_mt <- function(x)
{
  grepl('mt-',x)
}

directory <- "../upload/processeddata/"
sampleFiles_1 <- grep("hcg",list.files(directory),value=TRUE)
sampleName <- sub(".hcg.*","\\1",sampleFiles)
sampleCondition <- c("P815_M_dark","P815_M_dark","P815_M_dark","P815_M_light","P815_M_light",
                     "P815_M_light","P815_IFNgama_dark","P815_IFNgama_dark","P815_IFNgama_dark",
                     "P815_IFNgama_light","P815_IFNgama_light","P815_IFNgama_light")

# batch correct
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
keep <- !sapply(ddsHTSeq@rowRanges@partitioning@NAMES, is_mt)
ddsHTSeq <- ddsHTSeq[keep,]
mtx <- ddsHTSeq@assays@data$counts
batch <- c(1,1,1,1,1,1,1,1,0,1,1,0)
cov1 <- c(1,1,1,0,0,0,1,1,1,0,0,0)
cov2 <- c(0,0,0,0,0,0,1,1,1,1,1,1)
covar <- cbind(cov1,cov2)
adjusted_mtx <- ComBat_seq(mtx,batch = batch, covar_mod = covar)
mode(adjusted_mtx) <- "integer"


# f1 vs f2
index <- c(1,2,3,4,5,6)
sampleTable <- data.frame(sampleName = sampleName[index],
                          fileName = sampleFiles[index],
                          condition = sampleCondition[index])
sampleTable$condition <- factor(sampleTable$condition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
keep <- !sapply(ddsHTSeq@rowRanges@partitioning@NAMES, is_mt)
ddsHTSeq <- ddsHTSeq[keep,]
ddsHTSeq@assays@data$counts <- adjusted_mtx[,index]
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]
dds <- DESeq(ddsHTSeq)
res <- results(dds, filterFun=ihw)
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,
         col=colors,filename ='../result/f1_f2_sample.png')
select <- res$padj <= 0.05 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1
gene_f1_f2 <- rownames(res)[select]

# m2 vs m1
index <- c(7,8,9,10,11,12)
sampleTable <- data.frame(sampleName = sampleName[index],
                          fileName = sampleFiles[index],
                          condition = sampleCondition[index])
sampleTable$condition <- factor(sampleTable$condition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
keep <- !sapply(ddsHTSeq@rowRanges@partitioning@NAMES, is_mt)
ddsHTSeq <- ddsHTSeq[keep,]
ddsHTSeq@assays@data$counts <- adjusted_mtx[,index]
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]
dds <- DESeq(ddsHTSeq)
res <- results(dds, filterFun=ihw)
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,
         col=colors,filename ='../result/m1_m2_sample.png')
select <- res$padj <= 0.05 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1 &
  !(rownames(res) %in% gene_f1_f2)
gene_m1_m2 <- rownames(res)[select]
go_bp <- enrichGO(gene         = gene_m1_m2,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = 'SYMBOL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)
write.csv(go_bp@result, file="../result/m1_m2_go_bp.csv")
png(file='../result/m1_m2_go_bp.png',width=1200, height=1000,res=150)
dotplot(go_bp,showCategory = 10,x='count')
dev.off()
gene_go_bp <- c()
for (s in go_bp[1:10]$geneID) {
  g <- str_split(s,"/")[[1]]
  gene_go_bp <- c(gene_go_bp,g)
}
select_go_bp <- rownames(res) %in% gene_go_bp
df <- as.data.frame(colData(dds)[c("condition")])
pheatmap(assay(vsd)[select_go_bp,],cluster_rows = FALSE,
         color = colorRampPalette(c("blue","white","red"))(255),
         clustering_distance_cols = "correlation",annotation_col = df,
         filename ='../result/m1_m2_gene.png')

# f2 vs m2
index <- c(1,2,3,7,8,9)
sampleTable <- data.frame(sampleName = sampleName[index],
                          fileName = sampleFiles[index],
                          condition = sampleCondition[index])
sampleTable$condition <- factor(sampleTable$condition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
keep <- !sapply(ddsHTSeq@rowRanges@partitioning@NAMES, is_mt)
ddsHTSeq <- ddsHTSeq[keep,]
ddsHTSeq@assays@data$counts <- adjusted_mtx[,index]
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]
dds <- DESeq(ddsHTSeq)
res <- results(dds, filterFun=ihw)
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,
         col=colors,filename ='../result/f2_m2_sample.png')
select <- res$padj <= 0.05 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1
gene_f2_m2 <- rownames(res)[select]

# f1 vs m1
index <- c(4,5,6,10,11,12)
sampleTable <- data.frame(sampleName = sampleName[index],
                          fileName = sampleFiles[index],
                          condition = sampleCondition[index])
sampleTable$condition <- factor(sampleTable$condition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
keep <- !sapply(ddsHTSeq@rowRanges@partitioning@NAMES, is_mt)
ddsHTSeq <- ddsHTSeq[keep,]
ddsHTSeq@assays@data$counts <- adjusted_mtx[,index]
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]
dds <- DESeq(ddsHTSeq)
res <- results(dds, filterFun=ihw)
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,
         col=colors,filename ='../result/f1_m1_sample.png')
select <- res$padj <= 0.05 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1 &
  !(rownames(res) %in% gene_f2_m2)
gene_f1_m1 <- rownames(res)[select]
go_bp <- enrichGO(gene         = gene_f1_m1,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = 'SYMBOL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)
write.csv(go_bp@result, file="../result/f1_m1_go_bp.csv")
png(file='../result/f1_m1_go_bp.png',width=1200, height=1000,res=150)
dotplot(go_bp,showCategory = 10,x='count')
dev.off()
gene_go_bp <- c()
for (s in go_bp[1:10]$geneID) {
  g <- str_split(s,"/")[[1]]
  gene_go_bp <- c(gene_go_bp,g)
}
select_go_bp <- rownames(res) %in% gene_go_bp
df <- as.data.frame(colData(dds)[c("condition")])
pheatmap(assay(vsd)[select_go_bp,],cluster_rows = FALSE,
         color = colorRampPalette(c("blue","white","red"))(255),
         clustering_distance_cols = "correlation",annotation_col = df,
         filename ='../result/f1_m1_gene.png')