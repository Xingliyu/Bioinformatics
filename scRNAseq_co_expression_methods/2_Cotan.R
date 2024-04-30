library(Seurat)
library(dplyr)
library(patchwork)
library(data.table)
library(tidyverse)
library(COTAN)
library(Matrix)
library(ggrepel)
library(corrplot)
library(pheatmap)

mycolours <- c("A" = "#8491B4B2","B"="#E64B35FF")
my_theme = theme(axis.text.x = element_text(size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
                 axis.text.y = element_text( size = 14, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),
                 axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
                 axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"))

out_dir = "Data/"
panc_Seurat_Fib <- readRDS(file='panc_Seurat_Fib')
all.genes<-rownames(panc_Seurat_Fib)
panc <- panc_Seurat_Fib@assays$RNA@counts

count.data <- as.data.frame(as.matrix(panc))
obj = new("scCOTAN", raw = count.data)
obj = initRaw(obj, GEO = "GSE210347", sc.method = "10X", cond = "Fib")

#---------------------------------------------------------------------
t = "Fib_only"
print(paste("Condition ", t, sep = ""))

#---------------------------------------------------------------------
n_cells = length(colnames(obj@raw))
print(paste("n cells", n_cells, sep = " "))

if (!file.exists(out_dir)) {
  dir.create(file.path(out_dir))
}

if (!file.exists(paste(out_dir, "cleaning", 
                       sep = ""))) {
  dir.create(file.path(out_dir, "cleaning"))
}
n_it = 1

ttm = clean(obj)
obj = ttm$object

ttm$pca.cell.2

pdf(paste(out_dir, "cleaning/", t, "_", n_it, 
          "_plots_before_cells_exlusion.pdf", sep = ""))
ttm$pca.cell.2



ggplot(ttm$D, aes(x = n, y = means)) + geom_point() + 
  geom_text_repel(data = subset(ttm$D, 
                                n > (max(ttm$D$n) - 15)), aes(n, 
                                                              means, label = rownames(ttm$D[ttm$D$n > 
                                                                                              (max(ttm$D$n) - 15), ])), nudge_y = 0.05, 
                  nudge_x = 0.05, direction = "x", 
                  angle = 90, vjust = 0, segment.size = 0.2) + 
  ggtitle(label = "B cell group genes mean expression", 
          subtitle = " - B group NOT removed -") + 
  my_theme + theme(plot.title = element_text(color = "#3C5488FF", 
                                             size = 20, face = "italic", vjust = -10, 
                                             hjust = 0.02), plot.subtitle = element_text(color = "darkred", 
                                                                                         vjust = -15, hjust = 0.01))

dev.off()


nu_est = round(obj@nu, digits = 7)

plot.nu <- ggplot(ttm$pca_cells, aes(x = PC1, 
                                     y = PC2, colour = log(nu_est)))

plot.nu = plot.nu + geom_point(size = 1, 
                               alpha = 0.8) + scale_color_gradient2(low = "#E64B35B2", 
                                                                    mid = "#4DBBD5B2", high = "#3C5488B2", 
                                                                    midpoint = log(mean(nu_est)), name = "ln (nu)") + 
  ggtitle("Cells PCA coloured by cells efficiency") + 
  my_theme + theme(plot.title = element_text(color = "#3C5488FF", 
                                             size = 20), legend.title = element_text(color = "#3C5488FF", 
                                                                                     size = 14, face = "italic"), legend.text = element_text(color = "#3C5488FF", 
                                                                                                                                             size = 13), legend.key.width = unit(2, 
                                                                                                                                                                                 "mm"), legend.position = "right")

pdf(paste(out_dir, "cleaning/", t, "_plots_PCA_efficiency_colored.pdf", 
          sep = ""))
plot.nu
dev.off()
plot.nu

nu_df = data.frame(nu = sort(obj@nu), n = c(1:length(obj@nu)))

ggplot(nu_df, aes(x = n, y = nu)) + geom_point(colour = "#8491B4B2", 
                                               size = 1) + my_theme  #+ ylim(0,1) + xlim(0,70)

yset = 0.30  #threshold to remove low UDE cells
plot.ude <- ggplot(nu_df, aes(x = n, y = nu)) + 
  geom_point(colour = "#8491B4B2", size = 1) + 
  my_theme + ylim(0, 0.5) + xlim(0, 200) + 
  geom_hline(yintercept = yset, linetype = "dashed", 
             color = "darkred") + annotate(geom = "text", 
                                           x = 500, y = 0.5, label = paste("to remove cells with nu < ", 
                                                                           yset, sep = " "), color = "darkred", 
                                           size = 4.5)



pdf(paste(out_dir, "cleaning/", t, "_plots_efficiency.pdf", 
          sep = ""))
plot.ude

dev.off()
#> png 
#>   2

plot.ude


obj@meta[(nrow(obj@meta) + 1), 1:2] = c("Threshold low UDE cells:", 
                                        yset)

to_rem = rownames(nu_df[which(nu_df$nu < 
                                yset), ])

obj@raw = obj@raw[, !colnames(obj@raw) %in% 
                    to_rem]

ttm = clean(obj)

obj = ttm$object
ttm$pca.cell.2

nu_est = round(obj@nu, digits = 7)
plot.nu <- ggplot(ttm$pca_cells, aes(x = PC1, 
                                     y = PC2, colour = log(nu_est)))
plot.nu = plot.nu + geom_point(size = 2, 
                               alpha = 0.8) + scale_color_gradient2(low = "#E64B35B2", 
                                                                    mid = "#4DBBD5B2", high = "#3C5488B2", 
                                                                    midpoint = log(mean(nu_est)), name = "ln(nu)") + 
  ggtitle("Cells PCA coloured by cells efficiency: last") + 
  my_theme + theme(plot.title = element_text(color = "#3C5488FF", 
                                             size = 20), legend.title = element_text(color = "#3C5488FF", 
                                                                                     size = 14, face = "italic"), legend.text = element_text(color = "#3C5488FF", 
                                                                                                                                             size = 13), legend.key.width = unit(2, 
                                                                                                                                                                                 "mm"), legend.position = "right")

pdf(paste(out_dir, "cleaning/", t, "_plots_PCA_efficiency_colored_FINAL.pdf", 
          sep = ""))
plot.nu
dev.off()
#> png 
#>   2

plot.nu


obj@yes_yes = c()
obj = cotan_analysis(obj)

saveRDS(obj, file = paste(out_dir, t, ".cotan.RDS", sep = ""))

obj = get.coex(obj)

coex <- extract.coex(object = obj)

# saving the structure
saveRDS(obj, file = paste(out_dir, t, ".cotan.RDS", sep = ""))

saveRDS(coex, file = paste(out_dir, t, ".matrix.cotan.RDS", sep = ""))

obj = readRDS(paste(out_dir, t, ".cotan.RDS", sep = ""))

obj_coex = readRDS(paste(out_dir, t, ".matrix.cotan.RDS", sep = ""))

#==========================
# Just test               =
#==========================
t = "Fib_only"

Fibroblast_general = c("FAP", "LUM", "DCN", "PDPN", "COL1A1", "MFAP4")  #,'Hes3'
iCAFs = c("SOD2", "HAS1", "MT2A")
myCAFs = c("ACTA2", "TAGLN", "TPM2", "MYLK", "MYL9", "ACTG2", "POSTN", "SDC1", "MMP11", "CTHRC1")

p_value_Fib = get.pval(object = obj, gene.set.col = c(Fibroblast_general,iCAFs,myCAFs),
                         gene.set.row = c(Fibroblast_general,iCAFs,myCAFs))

coex <- extract.coex(object = obj)

COTAN_corr= coex[rownames(coex) %in% c(Fibroblast_general,iCAFs,myCAFs), 
         colnames(coex) %in% c(Fibroblast_general,iCAFs,myCAFs)]

pheatmap(COTAN_corr, cluster_rows=FALSE, cluster_cols=FALSE)

all_marker_genes <- c("FAP", "LUM", "DCN", "PDPN", "COL1A1", "MFAP4", "SOD2", "MT2A", "ACTA2", "TAGLN", "TPM2", "MYLK", "MYL9", "ACTG2", "POSTN", "SDC1", "MMP11", "CTHRC1")
COTAN_corr=COTAN_corr[all_marker_genes,all_marker_genes]

pheatmap(COTAN_corr, cluster_rows=FALSE, cluster_cols=FALSE)

#==========================
# Just test               =
#==========================

partial.coex.cotan = coex[rownames(coex) %in% c(Fibroblast_general,iCAFs,myCAFs),colnames(coex) %in% c(Fibroblast_general,iCAFs,myCAFs)]

#partial.coex.cotan = E16.5_cotan@coex[rownames(E16.5_cotan@coex) %in% c(tf1,tf2,hk),colnames(E16.5_cotan@coex) %in% c(tf1,tf2,hk)]


#partial.pval.cotan = p_value_E16.5[rownames(p_value_E16.5) %in% c(tf1,tf2,hk),colnames(p_value_E16.5) %in% c(tf1,tf2,hk)]
partial.pval.cotan = p_value_Fib
#partial.pval.cotan = partial.pval.cotan <= 0.05
#partial.coex.cotan[!partial.pval.cotan] <- 0

partial.coex.cotan = reshape2::melt(as.matrix(partial.coex.cotan))
colnames(partial.coex.cotan) = c("g1","g2","coex")
for (n in c(1:nrow(partial.coex.cotan))) {
  if (partial.coex.cotan[n,"g1"] == partial.coex.cotan[n,"g2"]) {
    partial.coex.cotan[n,"coex"]=0
  }
  
}



partial.coex.cotan$g1 <- factor(partial.coex.cotan$g1, c(Fibroblast_general,myCAFs,iCAFs))
partial.coex.cotan$g2 <- factor(partial.coex.cotan$g2, c(Fibroblast_general,myCAFs,iCAFs))

C = ggplot(partial.coex.cotan) + 
  geom_tile(aes(x=g1,y=g2, fill = coex),colour = "black", show.legend = TRUE) +
  #  facet_grid( g1 ~ g2  ,scales = "free", space = "free") + 
  scale_fill_gradient2(mid = "white",limits=c(round(min(partial.coex.cotan$coex),digits = 0), round(max(partial.coex.cotan$coex),digits = 0)),low = "#DC0000B2", high = "#3C5488B2")+
  #scale_fill_gradient2(low = "darkred", mid = "white",  high = "darkblue", midpoint = 0,na.value = "grey80", space = "Lab", guide = "colourbar", aesthetics = "fill", limits = lim_coex, oob=scales::squish)+ theme(legend.position="bottom")+
  theme(#legend.title = element_blank(),
    #strip.text.x = element_text(color = "red"),
    #axis.text.y = element_text(color = ),
    axis.text.x = element_text(angle=45,hjust=1,vjust=1.0),
    legend.position="bottom"
  )#+geom_text(aes(label=ifelse(t_hk == "hk", "H","")), color="grey", size=3)
C

#figure <- ggarrange(C, S,
#                   labels = c("Co.", "Sp."),
#                   ncol = 2, nrow = 1)
#figure

Fibroblast_general = c("FAP", "LUM", "DCN", "PDPN", "COL1A1", "MFAP4")  #,'Hes3'
iCAFs = c("SOD2", "HAS1", "MT2A")
myCAFs = c("ACTA2", "TAGLN", "TPM2", "MYLK", "MYL9", "ACTG2", "POSTN", "SDC1", "MMP11", "CTHRC1")

list.genes = list("myCAFs" = myCAFs,"Fibroblast_general" = Fibroblast_general,"iCAFs" = iCAFs)

plot_heatmap(df_genes = list.genes, sets = (1:3), conditions = "Fib_only", dir = "Data/")

#==========================
# Just test               =
#==========================
obj@meta

#--------------------------------------------------------------
ncol(obj@raw)

#--------------------------------------------------------------
as.numeric(obj@meta[3, 2]) - ncol(obj@raw)

#--------------------------------------------------------------
quant.p = get.GDI(obj)

#--------------------------------------------------------------
head(quant.p)

#--------------------------------------------------------------
Fibroblast_general = c("FAP", "LUM", "DCN", "PDPN", "COL1A1", "MFAP4")  #,'Hes3'
iCAFs = c("IL6", "IL11", "LIF", "CXCL1", "CXCL8", "SOD2", "HAS1", "MT2A")
myCAFs = c("ACTA2", "TAGLN", "TPM2", "MYLK", "MYL9", "MYH11", "ACTG2"
           , "POSTN", "SDC1", "MMP11", "CTHRC1")

#NPGs = c("Nes", "Vim", "Sox2", "Sox1", "Notch1", "Hes1", "Hes5", "Pax6")  #,'Hes3'
#PNGs = c("Map2", "Tubb3", "Neurod1", "Nefm", "Nefl", "Dcx", "Tbr1")
#hk = c("Calm1", "Cox6b1", "Ppia", "Rpl18", "Cox7c", "Erh", "H3f3a", "Taf1", "Taf2", 
#       "Gapdh", "Actb", "Golph3", "Mtmr12", "Zfr", "Sub1", "Tars", "Amacr")

text.size = 12



quant.p$highlight = with(quant.p, ifelse(rownames(quant.p) %in% 
                                           Fibroblast_general, "Fibroblast_general", ifelse(rownames(quant.p) %in%
                                                          iCAFs, "Constitutive", ifelse(rownames(quant.p) %in% 
                                                                              myCAFs, "myCAFs", "normal"))))

textdf <- quant.p[rownames(quant.p) %in% 
                    c(Fibroblast_general, myCAFs, iCAFs), ]

mycolours <- c(Constitutive = "#00A087FF", 
               Fibroblast_general = "#E64B35FF", myCAFs = "#F39B7FFF", 
               normal = "#8491B4B2")
f1 = ggplot(subset(quant.p, highlight == 
                     "normal"), aes(x = sum.raw.norm, y = GDI)) + 
  geom_point(alpha = 0.1, color = "#8491B4B2", 
             size = 2.5)
GDI_plot = f1 + geom_point(data = subset(quant.p, 
                                         highlight != "normal"), aes(x = sum.raw.norm, 
                                                                     y = GDI, colour = highlight), size = 2.5, 
                           alpha = 0.8) + geom_hline(yintercept = quantile(quant.p$GDI)[4], 
                                                     linetype = "dashed", color = "darkblue") + 
  geom_hline(yintercept = quantile(quant.p$GDI)[3], 
             linetype = "dashed", color = "darkblue") + 
  geom_hline(yintercept = 1.5, linetype = "dotted", 
             color = "red", size = 0.5) + scale_color_manual("Status", 
                                                             values = mycolours) + scale_fill_manual("Status", 
                                                                                                     values = mycolours) + xlab("log normalized counts") + 
  ylab("GDI") + geom_label_repel(data = textdf, 
                                 aes(x = sum.raw.norm, y = GDI, label = rownames(textdf), 
                                     fill = highlight), label.size = NA, 
                                 alpha = 0.5, direction = "both", na.rm = TRUE, 
                                 seed = 1234) + geom_label_repel(data = textdf, 
                                                                 aes(x = sum.raw.norm, y = GDI, label = rownames(textdf)), 
                                                                 label.size = NA, segment.color = "black", 
                                                                 segment.size = 0.5, direction = "both", 
                                                                 alpha = 0.8, na.rm = TRUE, fill = NA, 
                                                                 seed = 1234) + theme(axis.text.x = element_text(size = text.size, 
                                                                                                                 angle = 0, hjust = 0.5, vjust = 0.5, 
                                                                                                                 face = "plain", colour = "#3C5488FF"), 
                                                                                      axis.text.y = element_text(size = text.size, 
                                                                                                                 angle = 0, hjust = 0, vjust = 0.5, 
                                                                                                                 face = "plain", colour = "#3C5488FF"), 
                                                                                      axis.title.x = element_text(size = text.size, 
                                                                                                                  angle = 0, hjust = 0.5, vjust = 0, 
                                                                                                                  face = "plain", colour = "#3C5488FF"), 
                                                                                      axis.title.y = element_text(size = text.size, 
                                                                                                                  angle = 90, hjust = 0.5, vjust = 0.5, 
                                                                                                                  face = "plain", colour = "#3C5488FF"), 
                                                                                      legend.title = element_blank(), legend.text = element_text(color = "#3C5488FF", 
                                                                                                                                                 face = "italic"), legend.position = "bottom")  # titl)
legend <- cowplot::get_legend(GDI_plot)
GDI_plot = GDI_plot + theme(legend.position = "none")
GDI_plot



#--------------------------------------------------------------
Fibroblast_general = c("FAP", "LUM", "DCN", "PDPN", "COL1A1", "MFAP4")  #,'Hes3'
iCAFs = c("IL6", "IL11", "LIF", "CXCL1", "CXCL8", "SOD2", "HAS1", "MT2A")
myCAFs = c("ACTA2", "TAGLN", "TPM2", "MYLK", "MYL9", "MYH11", "ACTG2"
           , "POSTN", "SDC1", "MMP11", "CTHRC1")

#NPGs = c("Nes", "Vim", "Sox2", "Sox1", "Notch1", "Hes1", "Hes5", "Pax6")  #,'Hes3'
#PNGs = c("Map2", "Tubb3", "Neurod1", "Nefm", "Nefl", "Dcx", "Tbr1")
#hk = c("Calm1", "Cox6b1", "Ppia", "Rpl18", "Cox7c", "Erh", "H3f3a", "Taf1", "Taf2", 
#       "Gapdh", "Actb", "Golph3", "Mtmr12", "Zfr", "Sub1", "Tars", "Amacr")


list.genes = list("myCAFs" = myCAFs, "Fibroblast_general" = Fibroblast_general, 
                  "iCAFs" = iCAFs)

plot_heatmap(df_genes = list.genes, sets = c(1:3), 
             conditions = "Fib_only", dir = "Data/")


#--------------------------------------------------------------
print(sessionInfo())

#--------------------------------------------------------------
