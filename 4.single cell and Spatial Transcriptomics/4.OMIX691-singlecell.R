{
library(tidyverse) 
library(reshape2) 
library(ggsci) 
library(readxl) 
library(showtext) 
showtext_auto(enable=T) 
library(tibble) 
library(Hmisc) 
library(scales)
library(xlsx)
library(patchwork)
library(cowplot) 
library(Seurat) 
}
options(Seurat.object.assay.version = "v5") 

## 1.Preprocessing of Single-Cell Data ----
### 1.1.General Information ----
GSE = "OMIX691"
Info = "OMIX691-Vitiligo-scRNA-Homo"
# Code = "R codes were collected by Dr. Hu from the Second Xiangya Hospital of Central South University"

meta = read_xlsx( "OMIX691.xlsx", sheet = 1 )
colnames( meta )
meta = meta %>% 
  mutate( sample = gsub( "^.*ple ", "", `File Title` ) ) %>% 
  mutate( group = ifelse( grepl( "healthy", `File Title` ), "Healthy", "Vitiligo" ) )
head( meta )
save( GSE, Info, meta, file = "0-metadata.rdata" )

raw_path = data.frame(
  path = list.files( "./0-RAWDATA", full.names = T )
)  %>% 
  setNames( "path" ) %>% 
  mutate( id = gsub( "^.*/", "", path ) )  %>% 
  mutate( id = gsub( "\\..*$", "", id ) ) 

raw_path$sample = meta$sample[ match( raw_path$id, meta$`File ID` ) ]

scdata_ls = lapply( c(1:nrow(raw_path) ), function(xx){ 
  dt_dir = raw_path[ xx, "path" ] 
  dt_id = raw_path[ xx, "id" ]
  sample = raw_path[ xx, "sample" ]
  
  untar( dt_dir, exdir = "./0-RAWDATA/" )
  
  xxx =  Read10X( paste( "./0-RAWDATA/", sample, "/outs/filtered_gene_bc_matrices/GRCh38",sep = "" ) ) 
  
  yyy = CreateSeuratObject( counts = xxx, 
                            project = sample, 
                            min.cells = 2, 
                            min.features = 100 
  )
  
} )


class(scdata_ls)
sc_OMIX691 = merge( scdata_ls[[1]], scdata_ls[2:nrow(raw_path)] )

save( sc_OMIX691, file = "1.1-scOMIX691_raw_counts.rda" )
rm( scdata_ls, raw_path )

### 1.2.Seurat object structure and operations----
GSE = "OMIX691"
Info = "OMIX691-Vitiligo-scRNA-Homo"
load( "0-metadata.rdata" ) 

sc_OMIX691$sample = meta$sample[ match( sc_OMIX691$orig.ident, meta$sample ) ]
sc_OMIX691$group = ifelse( grepl( "H", sc_OMIX691$sample ), "healthy", "vitiligo" )

Idents( sc_OMIX691 ) = "orig.ident" 

DefaultAssay( sc_OMIX691 ) = "RNA" 

save( sc_OMIX691, file = "1.2-scOMIX691_raw_counts.rda" )

### 1.3.Initial quality control----

sc_OMIX691$Percent.Mito = PercentageFeatureSet(sc_OMIX691, pattern = "^MT-" ) 
hist( sc_OMIX691$Percent.Mito )

sc_OMIX691$Percent.Ercc = PercentageFeatureSet(sc_OMIX691, pattern = "^ERCC" ) 
hist( sc_OMIX691$Percent.Ercc )

sc_OMIX691$Percent.Ribo = PercentageFeatureSet(sc_OMIX691, pattern = "^RP[SL]" ) 
hist( sc_OMIX691$Percent.Ribo )


save( sc_OMIX691, file = "1.3-scOMIX691_raw_counts.rda" ) 

gc()


p1 <- FeatureScatter(object = sc_OMIX691, raster = T,
                     shuffle = T, pt.size = 1,
                     feature1 = "nFeature_RNA", feature2 = "nCount_RNA")+
  scale_color_igv() +
  guides(color = guide_legend(override.aes = list(size = 4))); p1

p2 = FeatureScatter(object = sc_OMIX691, raster = T,
                    shuffle = T, pt.size = 1,
                    feature1 = "Percent.Mito", feature2 = "Percent.Ribo")+
  scale_color_igv() +
  guides(color = guide_legend(override.aes = list(size = 4))); p2

p = p1 + p2 + plot_layout( guides='collect' ) # Patchwork for combining plots
p

ggsave( p, filename = "1.3_Expression_Data_Quality_Visualization.pdf", width = 10, height = 4.5 )
rm( p, p1, p2 )

### 1.4.Quality control: Remove low-quality cells and potential confounding factors----

hist( sc_OMIX691$nCount_RNA, breaks = 100 )
hist( sc_OMIX691$nFeature_RNA, breaks = 100 )
hist( sc_OMIX691$Percent.Mito, breaks = 100 )
hist( sc_OMIX691$Percent.Ribo, breaks = 100 )

table(sc_OMIX691$nCount_RNA > 2000 ) 
table(sc_OMIX691$nCount_RNA < 20000 ) 
table(sc_OMIX691$nFeature_RNA > 1000 )
table(sc_OMIX691$nFeature_RNA < 3000 ) 
table( sc_OMIX691$Percent.Mito < 20 ) 
table( sc_OMIX691$Percent.Ercc < 5 )
table( sc_OMIX691$Percent.Ribo < 45 )
sc_OMIX691 = subset( sc_OMIX691, Percent.Mito < 20 )
table(sc_OMIX691$orig.ident)
save( sc_OMIX691, file = "1.4_scProcessed.rda" )
gc()


### 1.5.Normalization and scaling----

library(future)
availableCores()
nbrOfWorkers()

plan()
plan("multisession", workers = 3 )
options (future.globals.maxSize = 2000 * 1024^4)

sc_OMIX691 <- NormalizeData(object = sc_OMIX691, 
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000,
                            verbose = T,
                            block.size = 1000 )

gc()

sc_OMIX691 <- ScaleData(object = sc_OMIX691,
                        
                        # vars.to.regress = c( "nCount_RNA", "Percent.Mito", "Percent.Ribo" ), # Remove the effect of these variables
                        verbose = T,
                        block.size = 500 )




### 1.6.Preliminary selection of highly variable genes----
sc_OMIX691 <- FindVariableFeatures(object = sc_OMIX691, 
                                   
                                   nfeatures = 2000 
)


VariableFeaturePlot(sc_OMIX691,
                    cols = c("grey", "steelblue"),
) 


VariableFeaturePlot(sc_OMIX691,
                    cols = c("grey", "steelblue"),
) +
  scale_y_log10()




save( sc_OMIX691, file = "1.6_scProcessed_Normalized_Scaled.rda" )




hist(colSums(sc_OMIX691, slot = "counts" ),
     breaks = 100,
     main = "Raw counts" )

hist(colSums(sc_OMIX691, slot = "data" ),
     breaks = 100,
     main = "Total expression after normalization" ) 

hist(colSums(sc_OMIX691, slot = "scale.data" ),
     breaks = 100,
     main = "Total expression after scaling" )


### 1.7.Dimension reduction and preliminary clustering----
sc_OMIX691 <- RunPCA(object = sc_OMIX691, 
                     verbose = TRUE, 
                     ndims.print = 1:5, 
                     nfeatures.print = 10,
                     seed.use = 124
)
gc()

ElbowPlot( sc_OMIX691, ndims = 50) 


sc_OMIX691 <- RunTSNE(object = sc_OMIX691, 
                      seed.use = 124, 
                      dims = 1:30,
                      verbose = TRUE
)


sc_OMIX691 <- RunUMAP(object = sc_OMIX691, 
                      seed.use = 124,
                      dims = 1:30,
                      verbose = TRUE
)

DimPlot( sc_OMIX691, reduction = "tsne" )
DimPlot( sc_OMIX691, reduction = "umap" )


p1 = PCAPlot( sc_OMIX691, 
              dims = c(1, 2),group.by = 'group',raster=FALSE)+
  ggtitle("PCA")+
  theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_color_lancet(); p1

#
p2 = TSNEPlot(sc_OMIX691, 
              dims = c(1, 2), group.by = 'group',raster=FALSE)+
  ggtitle("TSNE")+
  theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_color_lancet(); p2


p3 = UMAPPlot( sc_OMIX691, dims = c(1, 2), group.by = 'group',raster=FALSE)+
  ggtitle("UMAP")+
  theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_color_lancet(); p3


p = p1 + p2 + p3 + plot_layout( guides='collect' )
p
ggsave( p, filename = "1.7_Dimensionality_Reduction_Visualization.pdf", width = 10, height = 3.3 )

### 1.8.Clustering----

sc_OMIX691 <- FindNeighbors(sc_OMIX691, dims = 1:30) 

sc_OMIX691 <- FindClusters(object = sc_OMIX691,
                           resolution = 0.1, 
                           random.seed = 124, 
                           verbose = T
)

table( sc_OMIX691@meta.data$seurat_clusters)
save( sc_OMIX691, file = "1.8-scOMIX691_raw_counts.rda" )

DimPlot(sc_OMIX691, reduction = "pca", group.by = c("orig.ident", "unintegrated_clusters" ))
DimPlot(sc_OMIX691, reduction = "tsne", group.by = c("orig.ident", "unintegrated_clusters" ))
DimPlot(sc_OMIX691, reduction = "umap", group.by = c("orig.ident", "unintegrated_clusters" ))

p1 = TSNEPlot( sc_OMIX691,
               label =T,
               group.by = 'seurat_clusters',raster=FALSE)+
  ggtitle("TSNE")+
  theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_color_igv(); p1

#
p2 = UMAPPlot( sc_OMIX691,
               label =T,
               group.by = 'seurat_clusters',raster=FALSE)+
  ggtitle("UMAP")+
  theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_color_igv(); p2

#
p3 = UMAPPlot( sc_OMIX691,
               label =T,
               group.by = 'Group',raster=FALSE)+
  ggtitle("UMAP")+
  theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_color_igv(); p3
p1 + p2 + p3+ plot_layout( guides='collect' )
ggsave( filename = "1.8_Cell_Clustering_Visualization.pdf", width = 10, height = 4.8 )
rm(p, p1, p2 ,p3)


### 1.9.Marker gene analysis for clusters----
sc_OMIX691_int2 <- JoinLayers( sc_OMIX691 )
sc_OMIX691=sc_OMIX691_int2
sc.markers <- FindAllMarkers(object = sc_OMIX691, 
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             thresh.use = 0.5 
)

save( sc.markers, file = "1.9_AllMarkers.rda" )
library(openxlsx)
write.xlsx(sc.markers, file = "1.9_AllMarkers.xlsx" )


plotdata=sc_OMIX691@meta.data

top5 = sc.markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC)%>% 
  group_by()
write.xlsx(top5, file = "1.9_top5.xlsx" )

DoHeatmap(object = sc_OMIX691, 
          features = top5$gene, 
          label = TRUE
) +
  scale_fill_gradientn (colors = c ("steelblue", "white", "salmon")) 
ggsave( filename = "1.9_Cell_Clustering_and_Markers_Heatmap.pdf", width = 20, height = 15 )
dev.off()

DotPlot( object = sc_OMIX691, 
         group.by = "seurat_clusters",
         features = unique( top5$gene ) ) +
  coord_flip()

xxx = sc.markers %>% 
  filter( cluster == "7" )


ggsave( filename = "1.9_Cell_Clustering_and_Markers_Bubble_Plot.pdf", width = 10, height = 20 )

save( sc_OMIX691, file = "1.9_scProcessed_Merged_Clustering.rda" ) 

### 1.10.Cluster annotation----
plotdata=sc_OMIX691@meta.data

sc_OMIX691$cell_type = case_when( 
  sc_OMIX691$seurat_clusters %in% c( 4 ) ~ "Melanocyte",
  sc_OMIX691$seurat_clusters %in% c( 0,3,8,9,10,11,14 ) ~ "Keratinocyte",
  sc_OMIX691$seurat_clusters %in% c( 1 ) ~ "T cell",
  sc_OMIX691$seurat_clusters %in% c( 5) ~ "Langerhans cell",
  sc_OMIX691$seurat_clusters %in% c( 2,12) ~ "Mononuclear phagocyte",
  sc_OMIX691$seurat_clusters %in% c( 6 ) ~ "Fibroblast",
  sc_OMIX691$seurat_clusters %in% c( 7 ) ~ "Smooth muscle cell",
  sc_OMIX691$seurat_clusters %in% c( 13 ) ~ "Endothelial cell",
  sc_OMIX691$seurat_clusters %in% c( 15 ) ~ "B cell",
  sc_OMIX691$seurat_clusters %in% c( 16 ) ~ "Sebocytes",
)
save(sc_OMIX691,file="2.1_scProcessed_分群注释.rda")
## 2.single cell visualization code ----

### 2.1.Unannotated cluster plot ----
p1 = UMAPPlot( sc_OMIX691,
               label =F,
               group.by = 'seurat_clusters',raster=FALSE)+
  ggtitle("UMAP")+
  theme_bw()+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_fill_igv()+
  scale_color_d3('category20'); p1
pdf("2.1.Unannotated_cluster_plot.pdf",width=7,height=6)
print(p1)
dev.off()

### 2.2.Annotated cluster plot ----
p2 = UMAPPlot( sc_OMIX691,
               label =T,
               group.by = 'cell_type',raster=FALSE)+
  ggtitle("UMAP")+
  theme_bw()+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_fill_igv()+
  scale_color_d3('category20'); p2
pdf("2.2.Annotated_cluster_plot.pdf",width=7.5,height=6)
print(p2)
dev.off()

### 2.3.CTSS UMAP plot ----
p3 = FeaturePlot( sc_OMIX691, features = c("CTSS"), order = F,
                  split.by = "group", reduction = "umap")
pdf(file="2.3.CTSS_Expression_plot.pdf",width = 12,height = 10)
print(p1/p2|p3)
dev.off()

### 2.4.CTSS sample Violin plot ----
p4 = VlnPlot(melanocyte, features = c("CTSS"), group.by = "sample",split.plot = F, split.by = "Group",ncol=1
             ,stack = 2,raster = F,pt.size = 0,cols = c("#78bee5", "#e7a40e")) + 
  geom_boxplot(width = 0.3, outlier.shape = NA, position = position_dodge(width = 0.9),show.legend = F)+
  ggtitle("The expression levels of CTSS in melanocytes across different samples")+ 
  theme(plot.title = element_text(size = 12,face = "plain"))
p4
pdf(file="2.4.CTSS_expression_in_samples.pdf",width = 6,height=3)
print(p4)
dev.off()

### 2.5.CTSS proportion Violin plot ----
Idents(sc_OMIX691) <- "cell_type"
table(sc_OMIX691$orig.ident) 

prop.table(table(Idents(sc_OMIX691)))  # Proportion of cells in each cluster

table(Idents(sc_OMIX691), sc_OMIX691$sample)  # Count the cells for each cluster in each sample

Cellratio <- prop.table(table(Idents(sc_OMIX691), sc_OMIX691$sample), margin = 2)

Cellratio <- as.data.frame(Cellratio)
Cellratio$Var1 <- factor(Cellratio$Var1)
colourCount = length(unique(Cellratio$Var1))

# Bar plot showing the proportion of different cell types in each sample
library(ggplot2)
p5 = ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio',fill='')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  scale_fill_d3('category20')
pdf("Cell_type_proportion.pdf",width = 8,height = 6)
print(p5)
dev.off()

### 2.6.CTSS expression difference across different cell types ----
p6 = VlnPlot( 
  sc_OMIX691, features = c("CTSS"), 
  cols = c("#78bee5", "#e7a40e"),
  group.by = "cell_type",
  split.plot = F,
  split.by = "Group",
  ncol = 1,pt.size = 0)
p6
pdf("2.6.CTSS_expression_in_different_cell_types.pdf",width =7,height =4)
print(p6)
dev.off()

library(grid)
### 2.7.CTSS expression difference in melanocytes ----
p7 = VlnPlot( melanocyte, 
              features = c( "CTSS"), 
              cols = c("#78bee5", "#e7a40e"),
              flip = T, 
              fill.by = "ident" ,
              group.by = "Group",
              pt.size = 0
) +
  NoLegend()+
  ggtitle("Expression of CTSS in Melanocytes
        p = 2.94394E−48") +
  theme(plot.title = element_text(size = 10))
pdf("2.7.CTSS_expression_in_melanocytes.pdf",width =3.5,height =3.5)
print(p7)
dev.off()

### 2.8.Expression Ratio of CTSS+ and CTSS- Melanocytes ----
melanocyte <- subset(sc_OMIX691, cell_type == "Melanocyte")
DefaultAssay(melanocyte ) <- "RNA"
ctss_expression <- FetchData(melanocyte , vars = "CTSS")
melanocyte$CTSS_group <- ifelse(ctss_expression > 0, "CTSS+", "CTSS-")
table(melanocyte$CTSS_group)
Idents(melanocyte) <- "CTSS_group"
table(melanocyte$orig.ident)
prop.table(table(Idents(melanocyte)))
table(Idents(melanocyte), melanocyte$sample)
Cellratio <- prop.table(table(Idents(melanocyte), melanocyte$Group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio$Var1 <- factor(Cellratio$Var1, levels = c("CTSS+", "CTSS-"))
colourCount = length(unique(Cellratio$Var1))
p1 = ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='Sample', y = 'Ratio', fill='') +
  coord_flip() +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid")) +
  scale_fill_manual(values = c("CTSS+" = "#e7a40e", "CTSS-" = "#78bee5"))
Cellratio <- prop.table(table(Idents(melanocyte), melanocyte$sample), margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio$Var1 <- factor(Cellratio$Var1, levels = c("CTSS+", "CTSS-"))
colourCount = length(unique(Cellratio$Var1))
p2 = ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1), stat = "identity", width = 0.7, size = 0.5, colour = '#222222') + 
  theme_classic() +
  labs(x='Sample', y = 'Ratio', fill='') +
  coord_flip() +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid")) +
  scale_fill_manual(values = c("CTSS+" = "#e7a40e", "CTSS-" = "#78bee5"))

print(p1/p2)
ggsave(file="2.8_CTSS+_and_CTSS-_Melanocyte_Expression_Ratio_Plot.pdf", width=6, height = 6)
save(melanocyte, file="melanocyte.rda")

### 2.9.Melanocyte Subgroup Analysis in Disease Group (Vitiligo)----
{

  library(GEOmirror)
  library(GEOquery)
  library(affy)
  library(stringr)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(GSVA)
  library(GSEABase)
  library(DOSE)
  library(enrichplot)
  library(GSVA)
  library(ggplot2)
  library(ggrepel)
  library(ggsci)
  library(plyr)
  library(reshape2)
  library(dplyr)
  library(tidyverse)
  library(MetaboSignal) 
  library(viridis)
  library(RColorBrewer)
  library(circlize) 
  library(Cairo)
  library(stats)
  library(rstatix)
  library(xCell)
  library(dplyr)
  library(devtools)
  library(scales)
  library(gridExtra)
  library(dendextend)
  library(ggtree)
  library(ggfortify)
  library(cluster)
  library(limma)
  library(pheatmap)
  library(ggplotify)
  library(AnnoProbe)
  library(tidyr)
  library(tidyverse)
  library(e1071)
  library(parallel)
  library(HGNChelper)
  library(showtext) 
  library(readxl) 
  library(ggsignif) 
  library(tibble) 
  library(patchwork)
  showtext_auto(enable = T)
  
}
GSE="melanocyte_vitiligo"
melanocyte_vitiligo = subset(melanocyte, subset = group == "vitiligo")
p1 = TSNEPlot( melanocyte,
               label =T,
               group.by = 'CTSS_group',raster=FALSE)+
  ggtitle("TSNE")+
  theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_color_jama(); p1
degs_temp = FindMarkers( melanocyte_vitiligo, 
                             ident.1 = "High", ident.2 = "Low", 
                             group.by = "CTSS_group", 
                             min.pct = 0.2,
                             logfc.threshold = 0 # 因为得到的倍率都比较小，建议设置 0
)
degs_temp = degs_temp %>%
  mutate( DEG = factor( ifelse( abs(degs_temp$avg_log2FC) > 0.5 & degs_temp$p_val < 0.05,
                                ifelse(degs_temp$avg_log2FC > 0, 'up','down'),
                                'ns'), 
                        levels = c("up", "ns", "down" ), ordered = T ) ) %>% 
  mutate( symbol = rownames(.) )
table( degs_temp$DEG )
save(degs_temp, file = "vitiligo_MC_DEG.rds")
write.csv(degs_temp, file = "vitiligo_MC_DEG.csv")

### 2.10.Differential Expression Volcano Plot ----
logFC_cut=0.5
degs_temp=degs_temp[-1,]
p5 = ggplot(data = degs_temp,
            aes(x = avg_log2FC, y = -log10(p_val), color = DEG)) +
  geom_point(size = 1, shape = 16, alpha = 0.9) +
  theme_light() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_color_manual(values = c('#e7a40e','grey','#78bee5' ))+
  guides(color = guide_legend(override.aes = list(size = 4))) +
  geom_vline( xintercept = c(-logFC_cut, logFC_cut), linetype = "dashed" ) +
  geom_hline( yintercept = -log10(0.05), linetype = "dashed" )
p5
ggsave( p5, filename = paste("2.10._DEG火山图_",GSE,".pdf",sep= ""),width = 5, height = 4)
rm( logFC_cut )

degs_up = degs_temp %>%
  filter( DEG == "up" ) %>% 
  arrange( p_val) %>%
  pull( symbol ) %>% 
  bitr( ., fromType = "SYMBOL", toType = c("ENTREZID"),
        OrgDb = org.Hs.eg.db) 
degs_dn = degs_temp %>%
  filter( DEG == "down" ) %>% 
  arrange( p_val ) %>%
  pull( symbol ) %>%
  bitr( ., fromType = "SYMBOL", toType = c("ENTREZID"),
        OrgDb = org.Hs.eg.db) 

save( degs_up, degs_dn,file = paste("GSE","4.5_TopDEG.Rdata"))

### 2.11.TF Volcano Plot ----
tf_ls = read.table("/share/TF/Homo_sapiens_TF.txt", sep = "\t", header = T)
##xxx = checkGeneSymbols(tf_ls$Symbol)

de_tf = degs_temp %>%
  filter( symbol %in% tf_ls$Symbol ) %>% 
  mutate( symbol = ifelse( DEG == "ns", NA, symbol ) )
table( de_tf$DEG )

library(ggrepel)
p5 = ggplot(data = de_tf,
            aes(x = avg_log2FC, y = -log10(p_val), color = DEG)) +

  geom_point(size = 1, shape = 16, alpha = 0.9) +
  geom_text_repel( aes( label = symbol ), size = 3 ) +
  theme_light() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_color_manual(values = c('#e7a40e','grey','#78bee5' ))+
  guides(color = guide_legend(override.aes = list(size = 4))) +
  geom_vline( xintercept = c(-0.5, 0.5 ), linetype = "dashed" ) +
  geom_hline( yintercept = -log10(0.05), linetype = "dashed" )
p5
ggsave( p5, filename = paste("2.11_TF火山图_",GSE,".pdf",sep= ""),width = 5, height = 4)
rm( logFC_cut, p5 )


#
write.csv(de_tf, file = "2.11_差异TF.csv")



### 2.12.TF Prediction for CTSS ----
library(RcisTarget)
targene = list( targene = c(degs_up$SYMBOL, degs_dn$SYMBOL) )

data(motifAnnotations_hgnc)
motif <- importRankings("/share/TF/hg19-tss-centered-10kb-7species.mc9nr.feather")

tfpre = cisTarget( targene, motif, motifAnnot=motifAnnotations)

save( tfpre, file =  paste('2.12._转录因子预测_up',GSE,'.Rdata',sep = "") )

rm( tf_ls, de_tf, exr_tf, targene, motifAnnotations_hgnc, motif )

write.csv(tfpre,file="转录因子预测.csv")
tfpre=read.csv(file="转录因子预测.csv")

table( tfpre$TF_highConf == "" )
tfpre <- tfpre %>% 
  filter(grepl("CTSS", enrichedGenes))
#
plotdata = tfpre %>%
  filter( TF_highConf != "" ) %>% 
  mutate( TF_highConf = gsub("\\(.*\\). ","", TF_highConf ) ) %>% 
  mutate( TF = str_split(TF_highConf,"; ") ) %>%  
  unnest(TF) %>%  
  mutate( TF =  gsub(" ","", TF) ) %>% 
  mutate( Target = str_split(enrichedGenes,";") ) %>% 
  unnest( Target) %>% # 展开
  mutate( Target =  gsub(" ","", Target) ) %>% 
  distinct( TF,Target, .keep_all = T ) %>% 
  add_count( TF, name = "TF_c" ) %>% #
  add_count( Target, name = "Target_c" ) %>% 
  dplyr::select( TF, Target, NES, TF_c, Target_c ) 
data_filt = plotdata %>%
  filter( NES > quantile( NES, probs = 0 ) & 
            TF_c > quantile( TF_c, probs = 0 ) & 
            Target_c > quantile( Target_c, probs = 0 )
  ) 

ppi = data_filt %>% 
  dplyr::select( TF, Target, NES ) %>% 
  setNames( c("from", "to", "weight") )

nodes <- data.frame( degs = union( ppi$from, ppi$to ) )  %>% 
  mutate( Reg = degs_temp$avg_log2FC[ match( degs, degs_temp$symbol ) ] ) %>% 
  mutate( Reg = ifelse( Reg > 0 , "deg_up", "deg_down" ) ) %>% 
  mutate( Reg = ifelse( degs %in% ppi$from, "TF", Reg ) )

net <- igraph::graph_from_data_frame(d = ppi, vertices= nodes, directed = F)

igraph::V(net)$degree <- igraph::degree(net) 
# 
ggraph(net, layout = "fr", circular = FALSE) +
  geom_edge_fan(aes(edge_width = weight), color = "grey30", alpha = 0.3, show.legend = FALSE) +
  geom_node_point(aes(size = degree, color = Reg), alpha = 0.9) +
  geom_node_text(aes(size = degree, label = name), repel = TRUE) +  # 使用默认字体
  scale_edge_width(range = c(0.5, 1)) +
  scale_size_continuous(range = c(2, 10)) +
  theme_graph() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.background = element_rect(fill = "transparent", colour = NA),
    text = element_text(family = "")  # 移除所有特定字体设置，回到默认字体
  ) +
  scale_color_manual(values = c("#78bee5", "#e7a40e", "salmon3")) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(title = "Predicted TFs and Targets")

ggsave(filename = paste('2.12_预测TF的PPI_CTSS',GSE,'.pdf',sep = ""), width =10, height = 10 )

rm( plotdata, data_filt, ppi, nodes, net )

### 3.GO/KEGG enrichment
library(clusterProfiler)
library(stringr)
library(enrichplot)
library(data.table)
library(ggplot2)
IDs_up<-bitr(degs_up$SYMBOL,
             fromType='SYMBOL',
             toType=c('ENTREZID'),
             OrgDb='org.Hs.eg.db')

IDs_down<-bitr(degs_dn$SYMBOL,
               fromType='SYMBOL',
               toType=c('ENTREZID'),
               OrgDb='org.Hs.eg.db')
enrichGO_result_up<-enrichGO(IDs_up$ENTREZID,
                             OrgDb="org.Hs.eg.db",
                             qvalueCutoff=0.05,
                             pvalueCutoff=0.05,
                             ont="all")
enrichGO_result_down<-enrichGO(IDs_down$ENTREZID,
                               OrgDb="org.Hs.eg.db",
                               qvalueCutoff=0.05,
                               pvalueCutoff=0.05,
                               ont="all")

enrichGO_result_up = setReadable(enrichGO_result_up,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')
write.csv(enrichGO_result_up,file = "enrichGO_result_up.csv")
enrichGO_result_up=fread("enrichGO_result_up.csv")
enrichGO_result_down = setReadable(enrichGO_result_down,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')
write.csv(enrichGO_result_down,file = "eenrichGO_result_down.csv")
enrichGO_result_down=fread("eenrichGO_result_down.csv")
enrichGO_result_up_data<-data.frame(enrichGO_result_up)
enrichGO_result_down_data<-data.frame(enrichGO_result_down)

enrichGO_result_up_data<-mutate(enrichGO_result_up_data,up_down="Up-regulated")
enrichGO_result_up_data1=enrichGO_result_up_data[c(1:10),]

enrichGO_result_down_data<-enrichGO_result_down_data%>%
  mutate(up_down="Down-regulated",
         Count=-Count)
enrichGO_result_down_data1=enrichGO_result_down_data[c(1:5),]

enrichGO_result<-rbind(enrichGO_result_up_data1,enrichGO_result_down_data1)

label_position<-ifelse(enrichGO_result$up_down=="Down-regulated",0.5-enrichGO_result$Count,-0.5-enrichGO_result$Count)
label_hjust<-ifelse(enrichGO_result$up_down=="Down-regulated",0,1)


p<-ggplot(enrichGO_result,aes(x=Count,y=reorder(Description,Count),fill=up_down))+
  geom_bar(stat="identity",color="black",width=0.5)+
  geom_text(aes(label=str_wrap(Description,width=100)),position=position_nudge(x=label_position),hjust=label_hjust,size=3)+
  scale_fill_manual(values=c("#78bee5","#e7a40e"))+
  scale_x_continuous(limits=c(-30,30),breaks=seq(-30,30,by=10),labels=function(x)abs(x))+
  labs(x="Count",y="",title="")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.ticks.x=element_line(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        axis.title.x=element_text(size=12),
        legend.title=element_blank(),
        legend.position="top",#将图例位置设置在上面
        legend.text=element_text(size=12))
print(p)
ggsave("3.GO_barplot.pdf",width = 10,height = 8)

enrichKEGG_result_up<- enrichKEGG(gene = degs_up$ENTREZID, 
                                  keyType = "kegg",
                                  organism   = 'hsa',
                                  pvalueCutoff = 1,
                                  qvalueCutoff = 1,
                                  use_internal_data = F)
enrichKEGG_result_down<-enrichKEGG(gene = degs_dn$ENTREZID, 
                                   keyType = "kegg",
                                   organism   = 'hsa',
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 1,
                                   use_internal_data = F)
enrichKEGG_result_up = setReadable(enrichKEGG_result_up,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')
write.csv(enrichKEGG_result_up,file = "enrichKEGG_result_up.csv")

enrichKEGG_result_down = setReadable(enrichKEGG_result_down,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')
write.csv(enrichKEGG_result_down,file = "enrichKEGG_result_down.csv")
enrichKEGG_result_down=fread("enrichKEGG_result_down.csv")
enrichKEGG_result_up=fread("enrichKEGG_result_up.csv")

enrichKEGG_result_up<-data.frame(enrichKEGG_result_up)
enrichKEGG_result_down<-data.frame(enrichKEGG_result_down)



###加上一列：上调还是下调
enrichKEGG_result_up<-mutate(enrichKEGG_result_up,up_down="Up-regulated")
enrichKEGG_result_up1=enrichKEGG_result_up[c(1:5),]
enrichKEGG_result_down<-enrichKEGG_result_down%>%
  mutate(up_down="Down-regulated",
         Count=-Count)
enrichKEGG_result_down1=enrichKEGG_result_down
###两个数据框合并
enrichKEGG_result<-rbind(enrichKEGG_result_up1,enrichKEGG_result_down1)

###设置标签位置
label_position<-ifelse(enrichKEGG_result$up_down=="Down-regulated",0.5-enrichKEGG_result$Count,-0.5-enrichKEGG_result$Count)
label_hjust<-ifelse(enrichKEGG_result$up_down=="Down-regulated",0,1)


p<-ggplot(enrichKEGG_result,aes(x=Count,y=reorder(Description,Count),fill=up_down))+
  geom_bar(stat="identity",color="black",width=0.5)+
  geom_text(aes(label=str_wrap(Description,width=50)),position=position_nudge(x=label_position),hjust=label_hjust)+
  scale_fill_manual(values=c("#78bee5","#e7a40e"))+
  scale_x_continuous(limits=c(-30,30),breaks=seq(-30,30,by=10),labels=function(x)abs(x))+
  labs(x="Count",y="",title="")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.ticks.x=element_line(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        axis.title.x=element_text(size=12),
        legend.title=element_blank(),
        legend.position="top",
        legend.text=element_text(size=12))
print(p)
ggsave("3.0.kegg_barplot.pdf",width = 8,height = 6)

## 4.cellchat ----
{
  library(tidyverse)
  library(reshape2)
  library(ggsci)
  library(readxl) 
  library(showtext) 
  showtext_auto(enable=T)
  library(ggplotify) 
  library(tibble) 
  library(Hmisc)
  library(scales)
  library(xlsx)
  library(patchwork) 
  library(cowplot) 
  library(viridis)
  library(Seurat) # V 5.0.0
  library(CellChat) # V 1.6.1
  library( ggrepel )
  Info = "OMIX" 
}
head( sc_new@meta.data )
table( sc_new$group )

sc_new=sc_OMIX691

mtx_input = GetAssayData(sc_new, assay = "RNA", slot = "data")

mtx_c = mtx_input[ , sc_new$group == "healthy" ]
meta_c = sc_new@meta.data[ sc_new$group == "healthy", ]

mtx_e = mtx_input[ , sc_new$group == "vitiligo" ]
meta_e = sc_new@meta.data[ sc_new$group == "vitiligo", ]


cellchat_c <- createCellChat(object = mtx_c,
                             meta = meta_c,
                             group.by = "cell_type")

cellchat_e <- createCellChat(object = mtx_e, 
                             meta = meta_e,
                             group.by = "cell_type")

CellChatDB <- CellChatDB.human


cellchat_c@DB <- CellChatDB
cellchat_c <- subsetData(cellchat_c) 
plan()
future::plan("multicore", workers = 40 )
cellchat_c <- identifyOverExpressedGenes(cellchat_c)
cellchat_c <- identifyOverExpressedInteractions(cellchat_c)
cellchat_c <- projectData(cellchat_c, PPI.human)


cellchat_c <- computeCommunProb(cellchat_c)
cellchat_c <- aggregateNet(cellchat_c)

cellchat_c <- filterCommunication( cellchat_c, min.cells = 20)

cellchat_e@DB <- CellChatDB
cellchat_e <- subsetData(cellchat_e) 

cellchat_e <- identifyOverExpressedGenes(cellchat_e)
cellchat_e <- identifyOverExpressedInteractions(cellchat_e)
cellchat_e <- projectData(cellchat_e, PPI.human)


cellchat_e <- computeCommunProb(cellchat_e)
cellchat_e <- aggregateNet(cellchat_e)

cellchat_e <- filterCommunication( cellchat_e, min.cells = 20)


save( cellchat_c, cellchat_e, file = "4-分组Cellchat.rda" )
chat.list <- list( healthy = cellchat_c, vitiligo = cellchat_e)
cellchat_merge <- mergeCellChat( chat.list, add.names = names(chat.list), cell.prefix = TRUE )

save( cellchat_merge, file = "4-合并分组Cellchat.rda" )

pdf(file="4.1.正常组高低CTSS黑素细胞发出弦图.pdf",width = 12,height = 12)
netVisual_circle(cellchat_c@net$weight, 
                 edge.label.cex = 1.8,
                 
                 # edge.weight.max = 2, 
                 # edge.width.max = 2,
                 weight.scale = T, 
                 label.edge= T ,
                 vertex.label.cex = 1.4,
                 
                 arrow.size = 0.8,sources.use = c("MC_CTSS+","MC_CTSS-"))+title(main = "healthy sources")
dev.off()

pdf(file="4.2.疾病组高低CTSS黑素细胞发出弦图.pdf",width = 12,height = 12)
netVisual_circle(cellchat_e@net$weight, 
                 edge.label.cex = 1.8,
                 vertex.label.cex = 1.4,
                 # edge.weight.max = 2, 
                 # edge.width.max = 2,
                 weight.scale = T, 
                 label.edge= T,
                 arrow.size = 0.8,sources.use = c("MC_CTSS+","MC_CTSS-"))+title(main = "vitiligo sources")

dev.off()

pdf(file="4.1.正常组高低CTSS黑素细胞接收弦图.pdf",width = 12,height = 12)
netVisual_circle(cellchat_c@net$weight, 
                 edge.label.cex = 1.8,
                 vertex.label.cex = 1.4,
                 arrow.size = 0.8,
                 # edge.weight.max = 2, 
                 # edge.width.max = 2,
                 weight.scale = T,
                 
                 label.edge= T ,targets.use = c("MC_CTSS+","MC_CTSS-"))+title(main = "healthy targets")
dev.off()

pdf(file="4.2.疾病组高低CTSS黑素细胞接收弦图.pdf",width = 12,height = 12)
netVisual_circle(cellchat_e@net$weight, 
                 edge.label.cex = 1.8,
                 vertex.label.cex = 1.4,
                 arrow.size = 0.8,
                 # edge.weight.max = 2, 
                 # edge.width.max = 2,
                 weight.scale = T, 
                 label.edge= T,targets.use = c("MC_CTSS+","MC_CTSS-"))+title(main = "vitiligo targets")

dev.off()

pdf(file="4.3.疾病组chordMHC-I.pdf",width = 6,height = 6)
netVisual_chord_gene( cellchat_e, 
                      sources.use = c("MC_CTSS+","MC_CTSS-"), 
                      reduce = 0.005,
                      transparency = 0.8,
                      signaling = "MHC-I",
                      lab.cex = 0.8,
                      legend.pos.x = 6,
                      legend.pos.y = 6,
)+title("MHC-I vitiligo")


pdf(file="4.3.疾病组chordMHC-II.pdf",width = 6,height = 6)
netVisual_chord_gene( cellchat_e, 
                      sources.use = c("MC_CTSS+","MC_CTSS-"), 
                      reduce = 0.005,
                      transparency = 0.8,
                      signaling = "MHC-II",
                      lab.cex = 0.8,
                      legend.pos.x = 6,
                      legend.pos.y = 6,
                      
)+title("MHC-II vitiligo")
dev.off()


netVisual_bubble(cellchat_c, 
                 sources.use =c("MC_CTSS+","MC_CTSS-"),
                 targets.use = c("T cell","B cell","Mononuclear phagocyte","Langerhans cell"),
                 signaling = c('MHC-I',"MHC-II")  )
netVisual_bubble(cellchat_e, 
                 sources.use =c("MC_CTSS+","MC_CTSS-"),
                 targets.use = c("T cell","B cell","Mononuclear phagocyte","Langerhans cell"),
                 signaling = c('MHC-I',"MHC-II")  )

vitiligo=subset(sc_OMIX691,Group=="vitiligo")


vitiligo=subset(sc_OMIX691,Group=="vitiligo")
VlnPlot( vitiligo, 
         features = c( "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F","HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", 
                       "HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DMA","HLA-DMB" ), 
         stack = T, 
         flip = T, 
         fill.by = "ident" ,
         group.by = "cell_type"
) +
  ggtitle("MHC-I and MHC-II Expression in Vitiligo") +  # 设置图形标题
  theme(
    plot.title = element_text(hjust = 0.5),  # 标题水平居中
    strip.text = element_text(angle = 45),  # 调整标签文字角度
    strip.text.y = element_text(face = "plain", hjust = 0)  # 调整标签文字对齐方式
  ) +
  guides(fill = "none")
ggsave(file="4.4.疾病组小提琴MHC-I&MHC-II.pdf",width = 8,height = 10)



healthy=subset(sc_OMIX691,group=="healthy")
VlnPlot( healthy, 
         features = c( "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F","HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", 
                       "HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DMA","HLA-DMB" ), 
         stack = T, 
         flip = T, 
         fill.by = "ident" ,
         group.by = "cell_type"
) +
  ggtitle("MHC-I and MHC-II Expression in Healthy") +  # 设置图形标题
  theme(
    plot.title = element_text(hjust = 0.5),  # 标题水平居中
    strip.text = element_text(angle = 45),  # 调整标签文字角度
    strip.text.y = element_text(face = "plain", hjust = 0)  # 调整标签文字对齐方式
  ) +
  guides(fill = "none")
ggsave(file="4.4.正常组小提琴MHC-I&MHC-II.pdf",width = 8,height = 10)

VlnPlot( sc_OMIX691, 
         features = c( "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F","HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", 
                       "HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DMA","HLA-DMB" ), 
         split.plot = T,
         split.by = "Group",
         stack = T, 
         flip = T, 
         fill.by = "ident" ,
         group.by = "cell_type"
) +
  ggtitle("MHC-I and MHC-II Expression in Vitiligo") +  # 设置图形标题
  theme(
    plot.title = element_text(hjust = 0.5),  # 标题水平居中
    strip.text = element_text(angle = 45),  # 调整标签文字角度
    strip.text.y = element_text(face = "plain", hjust = 0)  # 调整标签文字对齐方式
  ) +
  guides(fill = "none")
ggsave(file="4.4.疾病组小提琴MHC-I&MHC-II.pdf",width = 8,height = 10)




pdf(file="4.气泡图.pdf",width = 12,height = 10)
netVisual_bubble(cellchat_merge,comparison = c(1,2),sources.use = c("MC_CTSS+","MC_CTSS-"),
                 targets.use = c("T cell","Mononuclear phagocyte","Langerhans cell"),
                 sort.by.target = TRUE)
dev.off()

netVisual_bubble(cellchat_merge,comparison = c(1,2),sources.use = c("MC_CTSS+","MC_CTSS-"),
                 targets.use = c("Keratinocyte"),
                 sort.by.target = TRUE)


pdf(file="4.疾病组KC到MC.pdf",width = 6,height = 6)
netVisual_chord_gene(cellchat_e, 
                     sources.use =c("Keratinocyte"),
                     targets.use = c("MC_CTSS+","MC_CTSS-"),
                     reduce = 0.005,
                     transparency = 0.8,
                     lab.cex = 0.8,
                     legend.pos.x = 6,
                     legend.pos.y = 6,
)+title("vitiligo")  
dev.off()
pdf(file="4.疾病组KC到MC气泡图.pdf",width = 6,height = 6)
netVisual_bubble(cellchat_merge,comparison = c(1,2),sources.use = c("Keratinocyte"),
                 targets.use = c("MC_CTSS+","MC_CTSS-"),
                 sort.by.target = TRUE)
dev.off()

## 5.ssGSEA analysis in single cell ----
library(GSVA)
melanocyte_vitiligo <- subset(melanocyte, subset = group == "vitiligo")
meta=melanocyte_vitiligo@meta.data
gene_data <- read_excel("16_cell_death.xlsx")
gene_lists <- list()

for (type in colnames(gene_data)) {
  genes <- as.character(gene_data[[type]])
  genes <- genes[!is.na(genes) & nchar(genes) > 0]
  gene_lists[[type]] <- genes
}

exp=as.matrix(GetAssayData(object = melanocyte_vitiligo, slot = "data"))
ssgsea <- gsva(exp, gene_lists, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
save(ssgsea,file = "ssgsea.rda")
a <- ssgsea %>% t() %>% as.data.frame()

identical(rownames(a), rownames(meta))
a$group <- meta$CTSS_group 
library(tibble)
a <- a %>% rownames_to_column("sample")
write.table(a, "ssGSEA.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
CTSS_ggsea <- gather(a,key = ssgsea, value = Expression, -c(group,sample)) 
p=ggplot(CTSS_ggsea, aes(x = ssgsea, y = Expression)) + 
  labs(y="Score", x =  NULL) +  
  geom_boxplot(aes(fill = group), position = position_dodge(0.5), width = 0.5, outlier.alpha = 0) + 
  scale_fill_manual(values = c("#78bee5", "#e7a40e")) +
  theme_bw() + 
  theme(plot.title = element_text(size = 12,color="black",hjust = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1,size=12 ),
        axis.text.y = element_text( size=12 ),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)) + 
  stat_compare_means(aes(group =  group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
p
ggsave(p,file="5.1-cell_death_differential_violin_plot.pdf",width=8,height = 7)


CTSS_exp=subset(exp,rownames(exp)=="CTSS")
CTSS_exp=t(CTSS_exp)
b= as.data.frame(CTSS_exp)
b <- b %>% rownames_to_column("sample")
a <- ssgsea %>% t() %>% as.data.frame()
a <- a %>% rownames_to_column("sample")
cor1=merge(b,a,by="sample")

rownames(cor1)=cor1[,1]
cor1=cor1[,-1]
cor1=t(cor1)

library(ggplot2)
library(reshape2)
cor1=t(cor1)
library(Hmisc)
cor_results <- rcorr(as.matrix(cor1), type="pearson")
cor_xx=cor_results$r
cor_CTSS=cor_xx[1,]
p=cor_results$P
p_CTSS=p[1,]
cor_melted <- melt(cor_CTSS)
colnames(cor_melted) <- c("Correlation")
cor_melted$Row="CTSS"
cor_melted$Column=rownames(cor_melted)
cor_melted=cor_melted[-1,]

p_melted <- melt(p_CTSS)
colnames(p_melted) <- c("PValue")
p_melted$Row="CTSS"
p_melted$Column=rownames(p_melted)
p_melted=p_melted[-1,]

data <- merge(cor_melted, p_melted)

data$Significance <- ifelse(data$PValue < 0.05, "*", "")
data$Significance <- ifelse(data$PValue < 0.0001, "****",
                            ifelse(data$PValue < 0.001, "***",
                                   ifelse(data$PValue < 0.01, "**",
                                          ifelse(data$PValue < 0.05, "*", ""))))
data$Label <- ifelse(data$Significance != "", sprintf("%.2f", data$Correlation), "")
p <- ggplot(data, aes(x=Column, y=Row, fill=Correlation)) +
  geom_tile() +  # 添加瓷砖图层
  geom_text(aes(label=Label), color="yellow", size=7, vjust=-0.5) +  # 显示相关系数，条件为显著时
  geom_text(aes(label=Significance), color="yellow", size=8, vjust=2) +  # 显示显著性标记，更大尺寸
  scale_fill_gradient2(low="#78bee5", high="brown", mid="white", midpoint=0, limit=c(-0, 0.3)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95,size = 18),
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y = element_text(size = 15),
        
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size=20,),
        legend.title = element_text(size = 15),  # 调整图例标题大小
        legend.text = element_text(size = 15)) + # 调整图例文本大小)+
  
  ggtitle("The correlations between 16 forms of cell death and CTSS in vitiligo melanocytes")
p
ggsave("5.2-CTSS_and_Cell_Death_Genes_Correlation_Heatmap.pdf", plot = p, width = 15, height = 4.5, units = "in") 



## 6.spatial transcriptomics ----
{
  library(tidyverse)
  library(reshape2)
  library(ggsci)
  library(readxl)
  library(showtext) 
  showtext_auto(enable=T)
  library(ggplotify) 
  library(tibble) 
  library(Hmisc)
  library(scales)
  library(xlsx)
  library(patchwork) 
  library(cowplot) 
  library(viridis)
  library(Seurat)
  library( ggrepel )
  library(xCell)
  library( GSVA )
  library(corrplot)
  library(Hmisc)
  Info = "Integrate"
  labxx = "Halo"
}
### 6.1.immune cells and cell death ssGSEA analysis ----
load("1.1_scHalo_Integrate.rda")
sct_matrix = as.matrix(scHalo@assays$SCT@scale.data)
sct_matrix[1:5,1:5]
xCell_result = xCellAnalysis( sct_matrix, rnaseq = F) 
xxx = CreateSeuratObject( xCell_result, project = "Immune", assay = "Immune" )
scHalo@assays$Immune = xxx@assays$Immune
save( scHalo, file = paste("1.1_scHalo_",Info, ".rda", sep = "" ) )
rm( xCell_result, xxx)
gene_data <- read_excel("./16_cell_death.xlsx")
gene_lists <- list()


for (type in colnames(gene_data)) {
  
  genes <- as.character(gene_data[[type]])
  genes <- genes[!is.na(genes) & nchar(genes) > 0]
  
  
  gene_lists[[type]] <- genes
}
sct_matrix = as.matrix(scHalo@assays$SCT@scale.data)
meta_result = gsva( sct_matrix, gene_lists, parallel.sz=20L, method =  "ssgsea" )

xxx = CreateSeuratObject( meta_result, project = "cell_death", assay = "cell_death" )
scHalo@assays$cell_death = xxx@assays$cell_death
expxx=as.matrix(GetAssayData(object = scHalo,assay = "cell_death",layer = "count"))
rownames(exp)
save(scHalo,file="1.1_scHalo_Integrate.rda")


### 6.2.痣中CTSS与细胞死亡相关性分析 ----
exp=as.matrix(GetAssayData(object = scHalo, assay = "cell_death",layer = "count"))
ctss=as.matrix(GetAssayData(object = scHalo, assay = "SCT",layer = "scale.data"))
identical(colnames(ctss),colnames(exp))
ctss=as.matrix(ctss["CTSS",])
colnames(ctss)="CTSS"
ctss=t(ctss)
expr=rbind(ctss,exp)
meta=scHalo@meta.data

meta_nevus=subset(meta,meta$Tissue=="Nevus")
nevus_expr = expr[, colnames(expr) %in% rownames(meta_nevus)]
nevus_expr=t(nevus_expr)

cor_results <- rcorr(as.matrix(nevus_expr), type="pearson")
cor_xx=cor_results$r
cor_CTSS=cor_xx[1,]
p=cor_results$P
p_CTSS=p[1,]

cor_melted <- melt(cor_CTSS)
colnames(cor_melted) <- c("Correlation")
cor_melted$Row="CTSS"
cor_melted$Column=rownames(cor_melted)
cor_melted=cor_melted[-1,]

p_melted <- melt(p_CTSS)
colnames(p_melted) <- c("PValue")
p_melted$Row="CTSS"
p_melted$Column=rownames(p_melted)
p_melted=p_melted[-1,]

data <- merge(cor_melted, p_melted)

data$Significance <- ifelse(data$PValue < 0.05, "*", "")




data$Significance <- ifelse(data$PValue < 0.0001, "****",
                            ifelse(data$PValue < 0.001, "***",
                                   ifelse(data$PValue < 0.01, "**",
                                          ifelse(data$PValue < 0.05, "*", ""))))


data$Label <- ifelse(data$Significance != "", sprintf("%.2f", data$Correlation), "")

p <- ggplot(data, aes(x=Column, y=Row, fill=Correlation)) +
  geom_tile() + 
  geom_text(aes(label=Label), color="yellow", size=7, vjust=-0.5) + 
  geom_text(aes(label=Significance), color="yellow", size=8, vjust=2) + 
  scale_fill_gradient2(low="#78bee5", high="brown", mid="white", midpoint=0, limit=c(-0.1, 0.62)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95,size = 18),
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y = element_text(size = 15),
        
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size=20,),
        legend.title = element_text(size = 15),  
        legend.text = element_text(size = 15)) +
  
  ggtitle("The correlations between 16 forms of cell death and CTSS in nevus")
p
ggsave("6.2-相关性热图.pdf", plot = p, width = 15, height = 4.5, units = "in") 


### 6.3.痣中CTSS与免疫细胞浸润相关性分析 ----
exp=as.matrix(GetAssayData(object = scHalo, assay = "Immune",layer = "count"))
ctss=as.matrix(GetAssayData(object = scHalo, assay = "SCT",layer = "scale.data"))
identical(colnames(ctss),colnames(exp))
ctss=as.matrix(ctss["CTSS",])
colnames(ctss)="CTSS"
ctss=t(ctss)
expr=rbind(ctss,exp)
meta=scHalo@meta.data
meta_nevus=subset(meta,meta$Tissue=="Nevus")
nevus_expr = expr[, colnames(expr) %in% rownames(meta_nevus)]
nevus_expr=t(nevus_expr)
nevus_expr=nevus_expr[,c("CTSS","aDC", "cDC","DC","Macrophages","Macrophages M1","Macrophages M2",
                         "B-cells", "NK cells","Tregs","CD4+ T-cells","CD8+ T-cells")]
cor_results <- rcorr(as.matrix(nevus_expr), type="pearson")
cor_xx=cor_results$r
cor_CTSS=cor_xx[1,]
p=cor_results$P
p_CTSS=p[1,]
cor_melted <- melt(cor_CTSS)
colnames(cor_melted) <- c("Correlation")
cor_melted$Row="CTSS"
cor_melted$Column=rownames(cor_melted)
cor_melted=cor_melted[-1,]
p_melted <- melt(p_CTSS)
colnames(p_melted) <- c("PValue")
p_melted$Row="CTSS"
p_melted$Column=rownames(p_melted)
p_melted=p_melted[-1,]
data <- merge(cor_melted, p_melted)

data$Significance <- ifelse(data$PValue < 0.05, "*", "")
data$Significance <- ifelse(data$PValue < 0.0001, "****",
                            ifelse(data$PValue < 0.001, "***",
                                   ifelse(data$PValue < 0.01, "**",
                                          ifelse(data$PValue < 0.05, "*", ""))))

data$Label <- ifelse(data$Significance != "", sprintf("%.2f", data$Correlation), "")

p <- ggplot(data, aes(x=Column, y=Row, fill=Correlation)) +
  geom_tile() + 
  geom_text(aes(label=Label), color="yellow", size=7, vjust=-0.5) + 
  geom_text(aes(label=Significance), color="yellow", size=8, vjust=2) +  
  scale_fill_gradient2(low="#78bee5", high="brown", mid="white", midpoint=0, limit=c(-0.2, 0.72)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95,size = 18),
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y = element_text(size = 15),
        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5,size=20,),
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15)) + 
  ggtitle("The correlations between  immune cell infiltration and CTSS in nevus") 
p
ggsave("6.3-免疫相关性热图.pdf", plot = p, width = 15, height = 3.4, units = "in") 

exp=as.matrix(GetAssayData(object = scHalo, assay = "Immune",layer = "count"))

ctss=as.matrix(GetAssayData(object = scHalo, assay = "SCT",layer = "scale.data"))
identical(colnames(ctss),colnames(exp))
ctss=as.matrix(ctss["CTSS",])
colnames(ctss)="CTSS"
ctss=t(ctss)
expr=rbind(ctss,exp)
meta=scHalo@meta.data

meta_nevus=subset(meta,meta$Tissue=="Nevus")
nevus_expr = expr[, colnames(expr) %in% rownames(meta_nevus)]
immune<- nevus_expr  %>%
  as.data.frame() %>%
  filter(rownames(.) %in% c("CTSS","aDC", "cDC","DC","Macrophages","Macrophages M1","Macrophages M2",
                            "B-cells", "NK cells","Tregs","CD4+ T-cells","CD8+ T-cells")) %>%
  t() %>%  
  as.data.frame() %>%  
  select(all_of(c("CTSS","aDC", "cDC","DC","Macrophages","Macrophages M1","Macrophages M2",
                  "B-cells", "NK cells","Tregs","CD4+ T-cells","CD8+ T-cells")))%>%
  select_if(~ sum(.) != 0) 
colnames(immune)=c("CTSS","aDC", "cDC","DC","Macrophages","M1","M2",
                   "B-cells", "NK cells","Tregs","CD4+ T","CD8+ T")
res <- rcorr(as.matrix(immune))  

cor_matrix <- res$r  
p_matrix <- res$P   

diag(p_matrix) <- 0
library(corrplot)
pdf("6.3-免疫浸润相关性热图.pdf",width = 11,height=11)
corrplot(cor_matrix, method = c('pie'), 
         type = c('upper'), 
         outline = 'grey', 
         
         diag = TRUE,
         tl.cex = 0.7, 
         tl.col = 'black',
         tl.pos = 'd',
         p.mat = p_matrix,
         col = colorRampPalette(c("dodgerblue4", "white", "brown"))(200),
         sig.level = c(.001, .01, .05),
         insig = "label_sig", 
         pch.cex = 1.2, 
         pch.col = 'black', 
         mar = c(0, 0, 0, 0) 
         
         
)

corrplot(cor_matrix, add = TRUE,
         method = c('number'), 
         type = c('lower'),
         
         
         diag = FALSE, 
         number.cex = 0.9,
         tl.pos = 'n', 
         cl.pos = 'n',
         col = colorRampPalette(c("dodgerblue4", "white", "brown"))(200),
         p.mat = p_matrix,
         insig = "pch",
         mar = c(0, 0, 0, 0) 
         
)
dev.off()

### 6.4.黑素指标，炎症指标，CTSS空间点图 ----
genexx = c( "TYR","DCT","MITF" )
tissue_chs = c( "Nevus" )

##
if(length(genexx)==1 ){
  namexx = "5.1-单指标-"
  plotdata = data.frame( scHalo@images$coord, 
                         Tissue = scHalo$Tissue,
                         Sample = scHalo$Sample,
                         Level = data_chs[ genexx, ]
  ) %>% 
    mutate( Tissue = ifelse( grepl( tissue_chs, Tissue ), 1, 0 )  ) 
  
  #
  ggplot( plotdata, aes( x = xt, y = yt, color = Level ) ) +
    geom_point( shape = 16, aes( size = Tissue ) ) +
    theme_light()+
    ggtitle( paste(genexx, sep = "") ) +
    theme(
      plot.title = element_text( hjust = 0.5),
      # legend.position = "bottom",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text( color = "black" ),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "transparent",colour = NA)
    )+
    scale_color_gradient( low = "gray90", high = "brown" ) +
    scale_size( range = c(1.2,2.2), breaks =  c(0,1), limits = c(0, 2), labels = c("Others", tissue_chs) )+
    xlim( cord_xx[tissue_chs, "x1"], cord_xx[tissue_chs,"x2"] )+
    ylim( cord_xx[tissue_chs,"y1"], cord_xx[tissue_chs,"y2"] )+
    facet_grid(Sample~.) 
  
}else{
  namexx = "5.2-多指标-"
  coordxx = scHalo@images$coord
  plotdata = lapply( genexx, function(xx){
    xxx = data_chs[ xx, ]
    xxx = data.frame( Expression = (xxx-min(xxx))/max(xxx-min(xxx)),
                      Gene = xx,
                      Tissue = scHalo$Tissue,
                      Sample = scHalo$Sample    ) %>% 
      mutate( Tissue = ifelse( grepl( tissue_chs, Tissue ), 1, 0 )  ) 
  }) %>% 
    do.call(rbind, . ) %>% 
    as.data.frame() %>% 
    mutate( yt = rep(coordxx$yt, length(genexx))  ) %>% 
    mutate( xt = rep(coordxx$xt, length(genexx))   ) %>% 
    group_by( yt ) %>% 
    mutate( Ratio = Expression/sum(Expression) ) %>% 
    mutate( Ratio = ifelse( Ratio == "NaN", 0, Ratio ) ) %>% 
    mutate( Gene = factor(Gene, levels = genexx, ordered = T) )
  
  #
  ggplot( plotdata, aes( x = xt, y = yt, color = Gene, alpha = Ratio*Expression ) ) +
    geom_point( shape = 16, aes( size = Tissue ) ) +
    theme_light()+
    ggtitle( "Levels" ) +
    theme(
      plot.title = element_text( hjust = 0.5),
      # legend.position = "bottom",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text( color = "black" ),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "transparent",colour = NA)
    )+
    scale_color_manual( values = rev( c("purple3","brown","#e7a40e","#78bee5","#e67c7d") ) ) +
    scale_alpha(range = c(0.05,1.2))+
    scale_size( range = c(1.2,2.2), breaks =  c(0,1), limits = c(0, 2), labels = c("Others", tissue_chs) )+
    xlim( cord_xx[tissue_chs,"x1"], cord_xx[tissue_chs,"x2"] )+
    ylim( cord_xx[tissue_chs,"y1"], cord_xx[tissue_chs,"y2"] )+
    facet_grid(Sample~.) +
    guides( alpha = "none" )
  
}

## 
ggsave( filename =paste( namexx, tissue_chs, "-", paste( sort(genexx), collapse = " & "), ".pdf", sep = "" ), 
        width =  cord_xx[tissue_chs, "wd"], height = cord_xx[tissue_chs, "ht"] ) # width = 5.9, 6.2




data_chs = do.call(rbind, 
                   list(
                     as.matrix( scHalo@assays[[ "SCT" ]]@scale.data ),
                     as.matrix( scHalo@assays[[ "Immune" ]]@data ),
                     as.matrix( scHalo@assays[[ "Function" ]]@data ),
                     as.matrix( scHalo@assays[[ "Metabolite" ]]@data ),
                     as.matrix( scHalo@assays[[ "Position" ]]@data )
                   ) )
#

cord_xx = data.frame( row.names = c( "Nevus",  "Peri nevus", "Epidermis", "Dermis", "Matrix" ),
                      x1 = c( 30, 35, 5, 2, 5),
                      x2 = c( 110, 120, 160, 157, 155),
                      y1 = c( 50, 37, 34, 30, 0 ),
                      y2 = c( 80, 80, 80, 76, 70),
                      wd = c( 4.2, 4.5, 5.8, 5.8, 5.4),
                      ht = c( 2.6, 3.5, 3, 3, 4.35 )
)

genexx = c( "CXCL9","CXCL10","CTSS" )
tissue_chs = c( "Nevus" )

##
if(length(genexx)==1 ){
  namexx = "5.1-单指标-"
  plotdata = data.frame( scHalo@images$coord, 
                         Tissue = scHalo$Tissue,
                         Sample = scHalo$Sample,
                         Level = data_chs[ genexx, ]
  ) %>% 
    mutate( Tissue = ifelse( grepl( tissue_chs, Tissue ), 1, 0 )  ) 
  
  #
  ggplot( plotdata, aes( x = xt, y = yt, color = Level ) ) +
    geom_point( shape = 16, aes( size = Tissue ) ) +
    theme_light()+
    ggtitle( paste(genexx, sep = "") ) +
    theme(
      plot.title = element_text( hjust = 0.5),
      # legend.position = "bottom",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text( color = "black" ),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "transparent",colour = NA)
    )+
    scale_color_gradient( low = "gray90", high = "brown" ) +
    scale_size( range = c(1.2,2.2), breaks =  c(0,1), limits = c(0, 2), labels = c("Others", tissue_chs) )+
    xlim( cord_xx[tissue_chs, "x1"], cord_xx[tissue_chs,"x2"] )+
    ylim( cord_xx[tissue_chs,"y1"], cord_xx[tissue_chs,"y2"] )+
    facet_grid(Sample~.) 
  
}else{
  namexx = "5.2-多指标-"
  coordxx = scHalo@images$coord
  plotdata = lapply( genexx, function(xx){
    xxx = data_chs[ xx, ]
    xxx = data.frame( Expression = (xxx-min(xxx))/max(xxx-min(xxx)),
                      Gene = xx,
                      Tissue = scHalo$Tissue,
                      Sample = scHalo$Sample    ) %>% 
      mutate( Tissue = ifelse( grepl( tissue_chs, Tissue ), 1, 0 )  ) 
  }) %>% 
    do.call(rbind, . ) %>% 
    as.data.frame() %>% 
    mutate( yt = rep(coordxx$yt, length(genexx))  ) %>% 
    mutate( xt = rep(coordxx$xt, length(genexx))   ) %>% 
    group_by( yt ) %>% 
    mutate( Ratio = Expression/sum(Expression) ) %>% 
    mutate( Ratio = ifelse( Ratio == "NaN", 0, Ratio ) ) %>% 
    mutate( Gene = factor(Gene, levels = genexx, ordered = T) )
  
  #
  ggplot( plotdata, aes( x = xt, y = yt, color = Gene, alpha = Ratio*Expression ) ) +
    geom_point( shape = 16, aes( size = Tissue ) ) +
    theme_light()+
    ggtitle( "Levels" ) +
    theme(
      plot.title = element_text( hjust = 0.5),
      # legend.position = "bottom",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text( color = "black" ),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "transparent",colour = NA)
    )+
    scale_color_manual( values = rev( c("purple3","brown","#e7a40e","#78bee5","#e67c7d") ) ) +
    scale_alpha(range = c(0.05,1.2))+
    scale_size( range = c(1.2,2.2), breaks =  c(0,1), limits = c(0, 2), labels = c("Others", tissue_chs) )+
    xlim( cord_xx[tissue_chs,"x1"], cord_xx[tissue_chs,"x2"] )+
    ylim( cord_xx[tissue_chs,"y1"], cord_xx[tissue_chs,"y2"] )+
    facet_grid(Sample~.) +
    guides( alpha = "none" )
  
}

## 
ggsave( filename =paste( namexx, tissue_chs, "-", paste( sort(genexx), collapse = " & "), ".pdf", sep = "" ), 
        width =  cord_xx[tissue_chs, "wd"], height = cord_xx[tissue_chs, "ht"] ) # width = 5.9, 6.2




data_chs = do.call(rbind, 
                   list(
                     as.matrix( scHalo@assays[[ "SCT" ]]@scale.data ),
                     as.matrix( scHalo@assays[[ "Immune" ]]@data ),
                     as.matrix( scHalo@assays[[ "Function" ]]@data ),
                     as.matrix( scHalo@assays[[ "Metabolite" ]]@data ),
                     as.matrix( scHalo@assays[[ "Position" ]]@data ),
                     as.matrix(GetAssayData(object = scHalo,assay = "cell_death",layer = "count"))
                   ) )
#
# unique( scHalo$Tissue )
cord_xx = data.frame( row.names = c( "Nevus",  "Peri nevus", "Epidermis", "Dermis", "Matrix" ),
                      x1 = c( 30, 35, 5, 2, 5),
                      x2 = c( 110, 120, 160, 157, 155),
                      y1 = c( 50, 37, 34, 30, 0 ),
                      y2 = c( 80, 80, 80, 76, 70),
                      wd = c( 4.2, 4.5, 5.8, 5.8, 5.4),
                      ht = c( 2.6, 3.5, 3, 3, 4.35 )
)

### 6.5.CTSS在单细胞白癜风和空间转录组与ICD相关基因的相关性分析 ----
exp=as.matrix(GetAssayData(object = scHalo, assay = "SCT",layer = "scale.data"))
meta=scHalo@meta.data
meta_nevus=subset(meta,meta$Tissue=="Nevus")
exp = exp[, colnames(exp) %in% rownames(meta_nevus)]
gene_data <- read_excel("./16_cell_death.xlsx")
gene_data<- gene_data %>%
  mutate(across(everything(), ~ sort(as.character(.), na.last = TRUE)))


valid_genes <- gene_data$`Immunogenic cell death`[gene_data$`Immunogenic cell death` %in% rownames(exp)]
Immunogenic_cell_death<- exp %>%
  as.data.frame() %>%
  filter(rownames(.) %in% c("CTSS", valid_genes)) %>%
  t() %>%  
  as.data.frame() %>%  
  select(all_of(c("CTSS", valid_genes)))%>%
  select_if(~ sum(.) != 0) 
res <- rcorr(as.matrix(Immunogenic_cell_death))  

cor_matrix <- res$r 
p_matrix <- res$P    
diag(p_matrix) <- 0

pdf("Immunogenic cell death.pdf", width = 6, height = 6)
corrplot(cor_matrix, 
         method = "pie",  
         type = "lower",  
         tl.srt = 45,     
         tl.col = "black",
         title = "Immunogenic cell death", 
         mar = c(0, 0, 1, 0), 
         col = colorRampPalette(c("dodgerblue4", "white", "brown"))(200), 
         p.mat = p_matrix,  
         sig.level = 0.05,  
         insig = "blank")  
dev.off()

load("./2.1_scProcessed_分群注释.rda")
xxx=subset(sc_OMIX691,Group == "vitiligo" & cell_type== "Melanocyte")
rm("sc_OMIX691")
exp1=as.matrix(GetAssayData(object = xxx, assay = "RNA",layer = "data"))
valid_genes1 <- gene_data$`Immunogenic cell death`[gene_data$`Immunogenic cell death` %in% rownames(exp1)]
valid_gene1=valid_genes
Immunogenic_cell_death1<- exp1 %>%
  as.data.frame() %>%
  filter(rownames(.) %in% c("CTSS", valid_genes)) %>%
  t() %>% 
  as.data.frame() %>%  
  select(all_of(c("CTSS", valid_genes)))%>%
  select_if(~ sum(.) != 0) 
res1 <- rcorr(as.matrix(Immunogenic_cell_death1)) 
cor_matrix1 <- res1$r  
p_matrix1 <- res1$P    
diag(p_matrix1) <- 0

pdf("6.5-免疫原性细胞死亡热图.pdf",width = 10,height=11)
corrplot(cor_matrix, method = c('pie'), 
         type = c('upper'), 
         outline = 'grey', 
         diag = TRUE,
         tl.cex = 1, 
         tl.col = 'black',
         tl.pos = 'd',
         p.mat = p_matrix,
         col = colorRampPalette(c("dodgerblue4", "white", "brown"))(200),
         sig.level = c(.001, .01, .05),
         insig = "label_sig", 
         pch.cex = 1.2, 
         pch.col = 'black', 
         mar = c(0, 0, 0, 0),
         cl.pos = "n"
)

corrplot(cor_matrix1, method = c('pie'), add = TRUE,
         type = c('lower'), 
         outline = 'grey', 
         diag = TRUE,
         tl.cex = 1.2, 
         tl.col = 'black',
         tl.pos = 'd',
         p.mat = p_matrix1,
         col = colorRampPalette(c("dodgerblue4", "white", "brown"))(200),
         sig.level = c(.001, .01, .05),
         insig = "label_sig",
         pch.cex = 1.2, 
         pch.col = 'black', 
         mar = c(0, 0, 0, 0) ,
)
mtext("The correlations between CTSS and immunogenic cell death genes expression in vitiligo melanocytes", side = 2, line = 1, las = 0, cex = 1.2)
mtext("The correlations between CTSS and immunogenic cell death genes expression in nevus", side = 3, line =0, las = 0, cex = 1.2)
dev.off()


### 6.6.空转荧光图 ----
data_chs = do.call(rbind, 
                   list(
                     as.matrix( scHalo@assays[[ "SCT" ]]@scale.data ),
                     as.matrix( scHalo@assays[[ "Immune" ]]@data ),
                     as.matrix( scHalo@assays[[ "Position" ]]@data ),
                     as.matrix(scHalo@assays$cell_death$counts)
                   ) )

cord_xx = data.frame( row.names = c( "Nevus",  "Peri nevus", "Epidermis", "Dermis", "Matrix" ),
                      x1 = c( 30, 35, 5, 2, 5),
                      x2 = c( 110, 120, 160, 157, 155),
                      y1 = c( 50, 37, 34, 30, 0 ),
                      y2 = c( 80, 80, 80, 76, 70),
                      wd = c( 4.2, 4.5, 5.8, 5.8, 5.4),
                      ht = c( 2.6, 3.5, 3, 3, 4.35 )
)

###########
genexx = c( "CTSS")
tissue_chs=c("Nevus")
dev.off()

{ 
 
  if(length(genexx)==1 ){
    namexx = "4.1-单指标-"
    plotdata = data.frame( scHalo@images$coord, 
                           Tissue = scHalo$Tissue,
                           Sample = scHalo$Sample,
                           Level = data_chs[ genexx, ]
    )
    plotdata <- plotdata %>%
      filter(xt > 20, xt <130, yt > 20, yt < 200)%>% 
      mutate( Tissue = ifelse( grepl( tissue_chs, Tissue ), 1, 0 )  ) 
    
    
    #
    ggplot( plotdata, aes( x = xt, y = yt, color = Level ) ) +
      geom_point( shape = 16, aes( size = Tissue ) ) +
      theme_light()+
      ggtitle( paste(genexx, sep = "") ) +
      theme(
        plot.title = element_text( hjust = 0.5),
        # legend.position = "bottom",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text( color = "black" ),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black",colour = NA),
        panel.background = element_rect(fill = "black", colour = NA),
        plot.margin = unit(c(0, 0, 0, 0), "lines"),
        legend.position = "none",  # 确认完全移除图例
        panel.border = element_blank()####分界线
      )+
      ylim(min(plotdata$yt) , max(plotdata$yt+1) )+
      xlim(min(plotdata$xt) , max(plotdata$xt) )+
      scale_size( range = c(2.3,2.3), breaks =  c(0,1), limits = c(0, 2), labels = c("Others", tissue_chs) )+
      scale_color_gradient( low = "black", high = "palegreen4",guide = F ) +
      coord_cartesian(expand = FALSE)+
      # scale_y_reverse()+
      
      facet_grid(Sample~.)
    
  }else{
    namexx = "4.2-多指标-"
    coordxx = scHalo@images$coord
    plotdata = lapply( genexx, function(xx){
      xxx = data_chs[ xx, ]
      xxx = data.frame( Expression = (xxx-min(xxx))/max(xxx-min(xxx)),
                        Gene = xx,
                        Sample = scHalo$Sample
      )  
    }) %>% 
      do.call(rbind, . ) %>% 
      as.data.frame() %>% 
      mutate( yt = rep(coordxx$yt, length(genexx))  ) %>% 
      mutate( xt = rep(coordxx$xt, length(genexx))   ) %>% 
      group_by( yt ) %>% 
      mutate( Ratio = Expression/sum(Expression) ) %>% 
      mutate( Ratio = ifelse( Ratio == "NaN", 0, Ratio ) ) %>% 
      mutate( Gene = factor(Gene, levels = genexx, ordered = T) )
    
    #
    ggplot( plotdata, aes( x = xt, y = yt, color = Gene, alpha = Ratio*Expression ) ) +
      geom_point( shape = 16, size = 1.5 ) +
      theme_light()+
      ggtitle( "Levels" ) +
      theme(
        plot.title = element_text( hjust = 0.5),
        # legend.position = "bottom",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text( color = "black" ),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)
      )+
      ylim( 0, 80 ) +
      facet_grid(Sample~.) +
      scale_alpha(range = c(0.01,1.2))+
      scale_color_manual( values = rev( c("orange3", "palegreen4","dodgerblue4","brown") ) ) +
      guides( alpha = "none" )
    
  }
} # 
## 
ggsave( filename =paste( namexx, paste( sort(genexx), collapse = " & "), ".tif", sep = "" ), 
        width = 5, height = 4.8 )

#########
genexx = c( "aDC")
tissue_chs=c("Nevus")
dev.off()
{ # 起始
  ##
  if(length(genexx)==1 ){
    namexx = "4.1-单指标-"
    plotdata = data.frame( scHalo@images$coord, 
                           Tissue = scHalo$Tissue,
                           Sample = scHalo$Sample,
                           Level = data_chs[ genexx, ]
    )
    plotdata <- plotdata %>%
      filter(xt > 20, xt <130, yt > 20, yt < 200)%>% 
      mutate( Tissue = ifelse( grepl( tissue_chs, Tissue ), 1, 0 )  ) 
    
    
    #
    ggplot( plotdata, aes( x = xt, y = yt, color = Level ) ) +
      geom_point( shape = 16, aes( size = Tissue ) ) +
      theme_light()+
      ggtitle( paste(genexx, sep = "") ) +
      theme(
        plot.title = element_text( hjust = 0.5),
        # legend.position = "bottom",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text( color = "black" ),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black",colour = NA),
        panel.background = element_rect(fill = "black", colour = NA),
        plot.margin = unit(c(0, 0, 0, 0), "lines"),
        legend.position = "none",  # 确认完全移除图例
        panel.border = element_blank()####分界线
      )+
      ylim(min(plotdata$yt) , max(plotdata$yt+1) )+
      xlim(min(plotdata$xt) , max(plotdata$xt) )+
      scale_size( range = c(2.3,2.3), breaks =  c(0,1), limits = c(0, 2), labels = c("Others", tissue_chs) )+
      scale_color_gradient( low = "black", high = "brown",guide = F ) +
      coord_cartesian(expand = FALSE)+
      # scale_y_reverse()+
      
      facet_grid(Sample~.)
    
  }else{
    namexx = "4.2-多指标-"
    coordxx = scHalo@images$coord
    plotdata = lapply( genexx, function(xx){
      xxx = data_chs[ xx, ]
      xxx = data.frame( Expression = (xxx-min(xxx))/max(xxx-min(xxx)),
                        Gene = xx,
                        Sample = scHalo$Sample
      )  
    }) %>% 
      do.call(rbind, . ) %>% 
      as.data.frame() %>% 
      mutate( yt = rep(coordxx$yt, length(genexx))  ) %>% 
      mutate( xt = rep(coordxx$xt, length(genexx))   ) %>% 
      group_by( yt ) %>% 
      mutate( Ratio = Expression/sum(Expression) ) %>% 
      mutate( Ratio = ifelse( Ratio == "NaN", 0, Ratio ) ) %>% 
      mutate( Gene = factor(Gene, levels = genexx, ordered = T) )
    
    #
    ggplot( plotdata, aes( x = xt, y = yt, color = Gene, alpha = Ratio*Expression ) ) +
      geom_point( shape = 16, size = 1.5 ) +
      theme_light()+
      ggtitle( "Levels" ) +
      theme(
        plot.title = element_text( hjust = 0.5),
        # legend.position = "bottom",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text( color = "black" ),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)
      )+
      ylim( 0, 80 ) +
      facet_grid(Sample~.) +
      scale_alpha(range = c(0.01,1.2))+
      scale_color_manual( values = rev( c("orange3", "palegreen4","dodgerblue4","brown") ) ) +
      guides( alpha = "none" )
    
  }
} # 
## 
ggsave( filename =paste( namexx, paste( sort(genexx), collapse = " & "), ".tif", sep = "" ), 
        width = 5, height = 4.8 )


##########
genexx = c( "CD4+ T-cells")
tissue_chs=c("Nevus")
dev.off()

{ # 起始
  ##
  if(length(genexx)==1 ){
    namexx = "4.1-单指标-"
    plotdata = data.frame( scHalo@images$coord, 
                           Tissue = scHalo$Tissue,
                           Sample = scHalo$Sample,
                           Level = data_chs[ genexx, ]
    )
    plotdata <- plotdata %>%
      filter(xt > 20, xt <130, yt > 20, yt < 200)%>% 
      mutate( Tissue = ifelse( grepl( tissue_chs, Tissue ), 1, 0 )  ) 
    
    
    #
    ggplot( plotdata, aes( x = xt, y = yt, color = Level ) ) +
      geom_point( shape = 16, aes( size = Tissue ) ) +
      theme_light()+
      ggtitle( paste(genexx, sep = "") ) +
      theme(
        plot.title = element_text( hjust = 0.5),
        # legend.position = "bottom",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text( color = "black" ),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black",colour = NA),
        panel.background = element_rect(fill = "black", colour = NA),
        plot.margin = unit(c(0, 0, 0, 0), "lines"),
        legend.position = "none",  # 确认完全移除图例
        panel.border = element_blank()####分界线
      )+
      ylim(min(plotdata$yt) , max(plotdata$yt+1) )+
      xlim(min(plotdata$xt) , max(plotdata$xt) )+
      scale_size( range = c(2.3,2.3), breaks =  c(0,1), limits = c(0, 2), labels = c("Others", tissue_chs) )+
      scale_color_gradient( low = "black", high = "brown",guide = F ) +
      coord_cartesian(expand = FALSE)+
      # scale_y_reverse()+
      
      facet_grid(Sample~.)
    
  }else{
    namexx = "4.2-多指标-"
    coordxx = scHalo@images$coord
    plotdata = lapply( genexx, function(xx){
      xxx = data_chs[ xx, ]
      xxx = data.frame( Expression = (xxx-min(xxx))/max(xxx-min(xxx)),
                        Gene = xx,
                        Sample = scHalo$Sample
      )  
    }) %>% 
      do.call(rbind, . ) %>% 
      as.data.frame() %>% 
      mutate( yt = rep(coordxx$yt, length(genexx))  ) %>% 
      mutate( xt = rep(coordxx$xt, length(genexx))   ) %>% 
      group_by( yt ) %>% 
      mutate( Ratio = Expression/sum(Expression) ) %>% 
      mutate( Ratio = ifelse( Ratio == "NaN", 0, Ratio ) ) %>% 
      mutate( Gene = factor(Gene, levels = genexx, ordered = T) )
    
    #
    ggplot( plotdata, aes( x = xt, y = yt, color = Gene, alpha = Ratio*Expression ) ) +
      geom_point( shape = 16, size = 1.5 ) +
      theme_light()+
      ggtitle( "Levels" ) +
      theme(
        plot.title = element_text( hjust = 0.5),
        # legend.position = "bottom",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text( color = "black" ),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)
      )+
      ylim( 0, 80 ) +
      facet_grid(Sample~.) +
      scale_alpha(range = c(0.01,1.2))+
      scale_color_manual( values = rev( c("orange3", "palegreen4","dodgerblue4","brown") ) ) +
      guides( alpha = "none" )
    
  }
} # 
## 
ggsave( filename =paste( namexx, paste( sort(genexx), collapse = " & "), ".tif", sep = "" ), 
        width = 5, height = 4.8 )

###############################
genexx = c( "CD8+ T-cells")
tissue_chs=c("Nevus")
dev.off()
{ # 起始
  ##
  if(length(genexx)==1 ){
    namexx = "4.1-单指标-"
    plotdata = data.frame( scHalo@images$coord, 
                           Tissue = scHalo$Tissue,
                           Sample = scHalo$Sample,
                           Level = data_chs[ genexx, ]
    )
    plotdata <- plotdata %>%
      filter(xt > 20, xt <130, yt > 20, yt < 200)%>% 
      mutate( Tissue = ifelse( grepl( tissue_chs, Tissue ), 1, 0 )  ) 
    
    
    #
    ggplot( plotdata, aes( x = xt, y = yt, color = Level ) ) +
      geom_point( shape = 16, aes( size = Tissue ) ) +
      theme_light()+
      ggtitle( paste(genexx, sep = "") ) +
      theme(
        plot.title = element_text( hjust = 0.5),
        # legend.position = "bottom",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text( color = "black" ),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black",colour = NA),
        panel.background = element_rect(fill = "black", colour = NA),
        plot.margin = unit(c(0, 0, 0, 0), "lines"),
        legend.position = "none",  # 确认完全移除图例
        panel.border = element_blank()####分界线
      )+
      ylim(min(plotdata$yt) , max(plotdata$yt+1) )+
      xlim(min(plotdata$xt) , max(plotdata$xt) )+
      scale_size( range = c(2.3,2.3), breaks =  c(0,1), limits = c(0, 2), labels = c("Others", tissue_chs) )+
      scale_color_gradient( low = "black", high = "brown",guide = F ) +
      coord_cartesian(expand = FALSE)+
      # scale_y_reverse()+
      
      facet_grid(Sample~.)
    
  }else{
    namexx = "4.2-多指标-"
    coordxx = scHalo@images$coord
    plotdata = lapply( genexx, function(xx){
      xxx = data_chs[ xx, ]
      xxx = data.frame( Expression = (xxx-min(xxx))/max(xxx-min(xxx)),
                        Gene = xx,
                        Sample = scHalo$Sample
      )  
    }) %>% 
      do.call(rbind, . ) %>% 
      as.data.frame() %>% 
      mutate( yt = rep(coordxx$yt, length(genexx))  ) %>% 
      mutate( xt = rep(coordxx$xt, length(genexx))   ) %>% 
      group_by( yt ) %>% 
      mutate( Ratio = Expression/sum(Expression) ) %>% 
      mutate( Ratio = ifelse( Ratio == "NaN", 0, Ratio ) ) %>% 
      mutate( Gene = factor(Gene, levels = genexx, ordered = T) )
    
    #
    ggplot( plotdata, aes( x = xt, y = yt, color = Gene, alpha = Ratio*Expression ) ) +
      geom_point( shape = 16, size = 1.5 ) +
      theme_light()+
      ggtitle( "Levels" ) +
      theme(
        plot.title = element_text( hjust = 0.5),
        # legend.position = "bottom",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text( color = "black" ),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)
      )+
      ylim( 0, 80 ) +
      facet_grid(Sample~.) +
      scale_alpha(range = c(0.01,1.2))+
      scale_color_manual( values = rev( c("orange3", "palegreen4","dodgerblue4","brown") ) ) +
      guides( alpha = "none" )
    
  }
} # 
## 
ggsave( filename =paste( namexx, paste( sort(genexx), collapse = " & "), ".tif", sep = "" ), 
        width = 5, height = 4.8 )