########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(cowplot)
library(kableExtra)
library(magrittr)
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
df_samples <- readxl::read_excel("doc/20191008_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% "test1")
df_samples <- df_samples[sample_n,]
print(df_samples)
(attach(df_samples))

#======1.2 load  SingleCellExperiment =========================
(load(file = "data/sce_4_20191008.Rda"))
names(sce_list)
object_list <- lapply(sce_list, as.Seurat)

for(i in 1:length(samples)){
    object_list[[i]]$orig.ident <- df_samples$c[i]
    object_list[[i]]$conditions <- df_samples$conditions[i]
    object_list[[i]]$group <- df_samples$group[i]
    object_list[[i]]$tests <- df_samples$tests[i]
}

#========1.3 merge ===================================
# https://support.bioconductor.org/p/76745/
# Don't use ComBat on raw counts
object <- Reduce(function(x, y) merge(x, y, do.normalize = F), object_list)
object@assays$RNA@data = object@assays$RNA@data *log(2)
remove(sce_list,object_list);GC()

(remove <- which(colnames(object@meta.data) %in% "ident"))
meta.data = object@meta.data[,-remove]
object@meta.data = meta.data 
remove(meta.data);GC()

#======== 1.4 FindVariableFeatures ===================================
# After removing unwanted cells from the dataset, the next step is to normalize the data.
#object <- NormalizeData(object = object, normalization.method = "LogNormalize", 
#                      scale.factor = 10000)
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(object), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
jpeg(paste0(path,"VariableFeaturePlot.jpeg"), units="in", width=10, height=7,res=600)
print(plot2)
dev.off()
#======1.5 1st run of pca-tsne  =========================
object <- ScaleData(object = object,features = VariableFeatures(object))
object <- RunPCA(object, features = VariableFeatures(object),verbose =F,npcs = 100)

npcs =100
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(reduction = "pca",resolution = 0.6,
                       dims.use = 1:npcs, print.output = FALSE)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs, check_duplicates = FALSE)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
object@assays$RNA@scale.data = matrix(0,0,0)

object@meta.data$orig.ident %<>% as.factor()
object@meta.data$orig.ident %<>% factor(levels = df_samples$sample)
p0 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 do.return = T,no.legend = F,label.size = 4, repel = T, title = "Original")
p1 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 no.legend = F,label.size = 4, repel = T, title = "Original")

#======1.6 Performing SCTransform and integration =========================
set.seed(100)
object_list <- SplitObject(object, split.by = "orig.ident")
remove(object);GC()
object_list %<>% lapply(SCTransform)
object.features <- SelectIntegrationFeatures(object_list, nfeatures = 3000)
options(future.globals.maxSize= object.size(object_list)*1.5)
object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = object.features, 
                                  verbose = FALSE)
anchors <- FindIntegrationAnchors(object_list, normalization.method = "SCT", 
                                  anchor.features = object.features)
object <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
object@assays$SCT@misc = list()
(load(file = "data/Integrated_KRPshS.RData"))
object@assays$integrated = combined@assays$integrated
object@reductions = combined@reductions
object@graphs = combined@graphs
object@neighbors = combined@neighbors
object@reductions$umap@cell.embeddings[,"UMAP_2"] = -object@reductions$umap@cell.embeddings[,"UMAP_2"]

remove(anchors,object_list);GC()
object %<>% RunPCA(npcs = 100, verbose = FALSE)
object <- JackStraw(object, num.replicate = 20,dims = 100)
object <- ScoreJackStraw(object, dims = 1:100)

jpeg(paste0(path,"JackStrawPlot.jpeg"), units="in", width=10, height=7,res=600)
JackStrawPlot(object, dims = 90:100)+
    ggtitle("JackStrawPlot")+
    theme(text = element_text(size=15),	
          plot.title = element_text(hjust = 0.5,size = 18))
dev.off()
npcs =20
object %<>% FindNeighbors(reduction = "pca",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.6)
object %<>% RunTSNE(reduction = "pca", dims = 1:npcs)
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
p2 <- TSNEPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 label.size = 4, repel = T,title = "Intergrated tSNE plot")
p3 <- UMAPPlot.1(object, group.by="orig.ident",pt.size = 1,label = F,legend.size = 15,
                 label.size = 4, repel = T,title = "Intergrated UMAP plot")
#=======1.9 summary =======================================
jpeg(paste0(path,"S1_remove_batch_tsne.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p0+ggtitle("Clustering without integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)),
          p2+ggtitle("Clustering with integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)))
dev.off()

jpeg(paste0(path,"S1_remove_batch_umap.jpeg"), units="in", width=10, height=7,res=600)
plot_grid(p1+ggtitle("Clustering without integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)),
          p3+ggtitle("Clustering with integration")+
              theme(plot.title = element_text(hjust = 0.5,size = 18)))
dev.off()
Idents(object) = "integrated_snn_res.0.6"
object %<>% sortIdent(numeric = T)
TSNEPlot.1(object = object, label = T,label.repel = T, group.by = "integrated_snn_res.0.6", 
         do.return = F, no.legend = F, title = "tSNE plot for all clusters, resolution = 0.6",
         pt.size = 0.3,alpha = 1, label.size = 5, do.print = T)

UMAPPlot.1(object = object, label = T,label.repel = T, group.by = "integrated_snn_res.0.5", 
           do.return = F, no.legend = F, title = "UMAP plot for all clusters, resolution = 0.5",
           pt.size = 0.2,alpha = 1, label.size = 5, do.print = T)

object@assays$integrated@scale.data = matrix(0,0,0)
save(object, file = "data/Coloratoral_3_20191203.Rda")
object_data = object@assays$RNA@data
save(object_data, file = "data/Coloratoral_data_4_20191008.Rda")

DefaultAssay(object) <- 'integrated'
for(i in c(4:16)/10){
    object %<>% FindClusters(resolution = i)
    Idents(object) = paste0("RNA_snn_res.",i)
    UMAPPlot.1(object, group.by=paste0("integrated_snn_res.",i),pt.size = 0.3,label = T,
               label.repel = T,alpha = 0.9,cols = Singler.colors,do.return = F,
               no.legend = F,label.size = 4, repel = T, title = paste("res =",i),
               do.print = T)
}
