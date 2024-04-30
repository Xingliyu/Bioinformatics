library(scLink)
library(caret)
library(corrplot)
library(pheatmap)
library(igraph)
library(corrr)
library(ggcorrplot)
library(lessR)

panc_Seurat_Fib <- readRDS(file='panc_Seurat_Fib')
panc <- t(as.matrix(panc_Seurat_Fib@assays$RNA@data))
class(panc)

#rownames(panc)
#colnames(panc)
panc[,1:3]
#genes[1:3,1:1]
#nchar(genes)

panc.norm = sclink_norm(panc, scale.factor = 1e6,  filter.genes = FALSE, gene.names = colnames(panc))
panc.norm = sclink_norm(panc, scale.factor = 1e6, filter.genes = TRUE, n = 1000)


any(is.na(panc.norm))

corr = sclink_cor(expr = panc.norm, ncores = 1)


#corrplot(corr, type = "upper", order = "hclust", 
#         tl.col = "black", tl.srt = 20,  col=colorRampPalette(c("blue","lightyellow","red"))(10))


#col <- colorRampPalette(c("blue", "white", "red"))
#corrplot(M, method="color", order="hclust", addrect=10, tl.pos="top", tl.cex = 0.5, tl.col = "black", col=col(100))

pheatmap(corr)
pheatmap(corr, cluster_cols = FALSE, cluster_rows = FALSE)

# sc_cor
#==========================================================================
corr.matrix0 <- corr
diag(corr.matrix0) <- 0

##set up threshold
threshold <- 0.7

## Now subsetting but here without absolute value
#ok <- apply(abs(corr.matrix0) >= threshold, 1, any)
ok <- apply( corr.matrix0 >= threshold, 1, any)

## or
# ok <- sort(unique( c(which(abs(corr.matrix0) >= threshold, arr = TRUE))))
# ok <- sort(unique( c(which(corr.matrix0 >= threshold, arr = TRUE))))

corr.matrixnew <-  corr[ok, ok]

pheatmap(corr.matrixnew)

Fibroblast_general = c("FAP", "LUM", "DCN", "PDPN", "COL1A1", "MFAP4")  #,'Hes3'
iCAFs = c("IL6", "IL11", "LIF", "CXCL1", "SOD2", "HAS1", "MT2A")
myCAFs = c("ACTA2", "TAGLN", "TPM2", "MYLK", "MYL9", "MYH11", "ACTG2"
           , "POSTN", "SDC1", "MMP11", "CTHRC1")

corr_maker_gene =  corr[rownames(corr) %in% c(Fibroblast_general,iCAFs,myCAFs), 
                        colnames(corr) %in% c(Fibroblast_general,iCAFs,myCAFs)]

pheatmap(corr_maker_gene, cluster_rows=FALSE, cluster_cols=FALSE)

all_marker_genes <- c("FAP", "LUM", "DCN", "PDPN", "COL1A1", "MFAP4", "SOD2", "MT2A", "ACTA2", "TAGLN", "TPM2", "MYLK", "MYL9", "ACTG2", "POSTN", "SDC1", "MMP11", "CTHRC1")
corr_maker_gene=corr_maker_gene[all_marker_genes,all_marker_genes]

pheatmap(corr_maker_gene, cluster_rows=FALSE, cluster_cols=FALSE)
#==========================================================================


# Make an Igraph object from this matrix:
network <- graph_from_adjacency_matrix(corr_maker_gene, weighted=T, mode="undirected", diag=F)

networks = sclink_net(expr = panc.norm, ncores = 1, lda = seq(0.5, 0.1, -0.05))


networks$cor[1:3,1:3]

# sc_net
#==========================================================================
corr.matrix0 <- networks$cor
diag(corr.matrix0) <- 0

##set up threshold
threshold <- 0.7
ok <- apply( corr.matrix0 >= threshold, 1, any)

corr.matrixnew <-  networks$cor[ok, ok]

pheatmap(corr.matrixnew)

Fibroblast_general = c("FAP", "LUM", "DCN", "PDPN", "COL1A1", "MFAP4")  #,'Hes3'
iCAFs = c("IL6", "IL11", "LIF", "CXCL1", "CXCL8", "SOD2", "HAS1", "MT2A")
myCAFs = c("ACTA2", "TAGLN", "TPM2", "MYLK", "MYL9", "MYH11", "ACTG2"
           , "POSTN", "SDC1", "MMP11", "CTHRC1")

net_corr = networks$cor

corr_maker_gene =  net_corr[rownames(net_corr) %in% c(Fibroblast_general,iCAFs,myCAFs), 
                        colnames(net_corr) %in% c(Fibroblast_general,iCAFs,myCAFs)]

pheatmap(corr_maker_gene)
