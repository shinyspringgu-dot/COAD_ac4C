setwd("~/mw/CRC")
sce.all = readRDS("sce.all.RDS")
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
sce.all <- subset(sce.all, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
sce.all <- NormalizeData(sce.all)
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst", nfeatures = 2000)
sce.all <- ScaleData(sce.all, features = rownames(sce.all))
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all),reduction.name = "pca")
ElbowPlot(sce.all)
library(harmony)
sce.all <- RunHarmony(sce.all,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")

sce.all <- FindNeighbors(sce.all, reduction = "harmony", dims = 1:15)
### 设置多个resolution选择合适的resolution
sce.all <- FindClusters(sce.all, resolution = 0.8)
levels(sce.all)
cols = c("#FB8072","#7BAFDE","#B17BA6","#FDB462","#E78AC3","#B2DF8A","#8DD3C7","#E6AB02","#BEAED4","#d4b7b7","#ba5ce3",
         "#aeae5c","#00bfff","#56ff0d", "#ffff00","#E41B1B", "#4376AC", "#48A75A", "#87638F", "#D87F32", "#737690", "#D690C6", "#B17A7D", "#847A74", "#4285BF",
         "#204B75", "#588257", "#B6DB7B", "#E3BC06", "#FA9B93", "#E9358B", "#A0094E", "#999999", "#6FCDDC", "#BD5E95")
sce.all@meta.data = sce.all@meta.data %>% 
  dplyr::mutate(celltype= case_when(
    seurat_clusters %in% c(7,15,17) ~ "Fibroblasts",
    seurat_clusters %in% c(13,20) ~ "Endothelials",
    seurat_clusters %in% c(6,18) ~ "B cells",
    seurat_clusters %in% c(3,23) ~ "Plasma",
    seurat_clusters %in% c(1,2,4,5,11) ~ "T cells",
    seurat_clusters %in% c(8,9) ~ "Myeloid",
    seurat_clusters %in% c(21) ~ "Mast cells",
    TRUE ~ "Epithelials"
  ))
DimPlot(sce.all,cols=cols,label=T,group.by="celltype")
DimPlot(sce.all,cols=cols,label=T)


####figure1C
Idents(sce.all) = sce.all$celltype
library(RColorBrewer) 
library(viridis)
library(wesanderson)
levels(sce.all)
#pal <- wes_palette("Zissou1", 10, type = "continuous")
pal = viridis(option="magma",direction = -1,n=15)
genes = c(#"CENPF","HMGB2",
  
  
  "EPCAM","KRT19",
  "LYZ","C1QA",
  "PLVAP","PECAM",
  "JCHAIN","IGKC",
  
  "CD3D","CD3E",
  "TPSB2","TPSAB1",
  "COL1A1","COL1A2",
  "MS4A1","CD79A"
  
  
  #"NKG7","KLRG1",
  
  
  
  #"TAGLN","ACTA2",
  #"GUCA2A","GUCA2B",
  
  
)
DotPlot(sce.all, features = genes ) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = c("#7BAFDE","white","#DC050C"))+ theme(legend.position="right")  + labs(title = "", y = "", x="")
