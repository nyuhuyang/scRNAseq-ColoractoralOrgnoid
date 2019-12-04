library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(fgsea)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
marker_path <- paste0(path,"markers/")
if(!dir.exists(marker_path))dir.create(marker_path, recursive = T)

# ======== 2.1 =========== test with known markers==================
(load(file = "data/Coloratoral_3_20191203.Rda"))
DefaultAssay(object) <- "RNA"
df_markers <- readxl::read_excel("../seurat_resources/bio-rad-markers.xlsx",sheet = "Mouse.Main")
colnames(df_markers) = gsub(" ","_",colnames(df_markers))
colnames(df_markers) = gsub(":|\\/","_",colnames(df_markers))
colnames(df_markers) = gsub("\\+","",colnames(df_markers))
markers = df_markers[,grep("Alias",colnames(df_markers),invert = T)]
marker.list <- df2list(markers)

marker.list %<>% lapply(function(x) x[1:18]) %>% 
     lapply(function(x) FilterGenes(object,x)) %>% 
     lapply(function(x) x[!is.na(x)]) %>% 
    lapply(function(x) x[1:min(length(x),12)])
marker.list <- marker.list[!is.na(marker.list)]
marker.list %>% list2df %>% t %>% kable() %>% kable_styling()

Idents(object) = "integrated_snn_res.0.6"
for(i in 1:length(marker.list)){
    if(length(marker.list[[i]]) == 0) next
    p <- lapply(marker.list[[i]], function(marker) {
        FeaturePlot(object = object, feature = marker,pt.size = 0.5,
                    reduction="umap", label = T)+
            NoLegend()+
            ggtitle(paste0(marker,Alias(df = df_markers,gene = marker)))+
            theme(plot.title = element_text(hjust = 0.5,size = 15,face = "plain"))
    })
    jpeg(paste0(path,"markers/",names(marker.list)[i],".jpeg"),units="in", width=10, height=7,res=600)
    print(do.call(cowplot::plot_grid, p)+ ggtitle(paste(names(marker.list)[i],"markers"))+
              theme(plot.title = element_text(hjust = 0.5,size = 20)))
    dev.off()
    Progress(i,length(marker.list))
}

## PLOT RELATIVE GENE EXPRESSION IN FEATUREPLOT

# Select multiple genes of interest
Endo_set <- c("Chga","Chgb","Syp")
Ent_set <- c("Krt20","Vil1","Slc26a3")
Fetal_set <- c("Krt4", "Spp1", "Ly6a")
Goblet_set <- c("Agr2","Spink4","Tff3")
Paneth_set <- c("Lyz1","Defa24","Defa2")
Reg_Pan <- c("Gm15284","AY761184","Defa17","Gm14851","Defa22","Defa-rs1",
             "Defa3","Defa24","Defa26","Defa21","Lyz1","Gm15292",
             "Mptx2","Ang4")
revSC_set <- c("Clu","Anxa1","Basp1")
Stem_set <- c("Lgr5","Smoc2","Olfm4")
Yap_set <- c("Ankrd1","Ccn1","Ccn2")

marker.list <- list("Endo_set" = Endo_set,
                    "Ent_set"=Ent_set,
                    "Fetal_set"=Fetal_set,
                    "Goblet_set" = Goblet_set,
                    "Paneth_set"=Paneth_set,
                    "Reg_Pan" = Reg_Pan,
                    "Stem_set"=Stem_set,
                    "revSC_set"=revSC_set,
                    "Yap_set" = Yap_set)

DefaultAssay(object) <- "SCT"
marker.list %<>%     lapply(function(x) FilterGenes(object,x)) 
marker.list %>% list2df %>% t %>% kable() %>% kable_styling()

for(i in seq_along(marker.list)){
    object[[names(marker.list)[i]]] <- PercentageFeatureSet(object,
                                                            features = marker.list[[i]])
    Progress(i, length(marker.list))
}

# Create custom color palette
GrYlOrRd <- c("lightgrey","#FFFFCC","#FFEDA0","#FED976","#FEB24C","#FD8D3C",
              "#FC4E2A","#E31A1C","#BD0026","#800026")

# Plot mean expression using Seurat::FeaturePlot()
object@meta.data$orig.ident %<>% as.factor()
object@meta.data$orig.ident %<>% factor(levels = c("PRE","POST","RES"))
table(object$orig.ident)
Idents(object) = "integrated_snn_res.0.6"

for(i in seq_along(marker.list)){
    FeaturePlot.1(object = object,
                  reduction = "umap",
                  cols = GrYlOrRd,
                  features = names(marker.list)[i], 
                  split.by = "orig.ident",
                  label = T,
                  pt.size = 0.5,
                  label.size = 4,
                  strip.text.size = 15,
                  title = names(marker.list)[i],
                  unique.name = T,
                  do.print = T) 
    Progress(i, length(marker.list))
}

#======== rename ident =================
Idents(object) = "integrated_snn_res.0.6"
object %<>% RenameIdents("0" = "Fetal cells",
                         "1" = "Enterocytes",
                         "2" = "Fetal cells",
                         "3" = "Stem cells",
                         "4" = "Stem cells",
                         "5" = "Stem cells",
                         "6" = "Fetal cells",
                         "7" = "Goblet cells",
                         "8" = "Unknown",
                         "9" = "Revival stem cells",
                         "10" = "Stem/Fetal cells",
                         "11" = "Unknown",
                         "12" = "Enteroendocrine")

object@meta.data$cell.types = as.character(Idents(object))
Idents(object) = "cell.types"
object %<>% sortIdent()
object <- AddMetaColor(object = object, label= "cell.types",
                       colors = c(Singler.colors,Singler.colors))
for( label in c(T,F)){
    UMAPPlot.1(object, group.by = "cell.types",cols = ExtractMetaColor(object),label = label,
               label.repel = T, pt.size = 0.5,label.size = 5, repel = T,no.legend = T,
               do.print = T,do.return = F,title = "Cell types in all 3 samples")
}
UMAPPlot.1(object, group.by = "cell.types",cols = ExtractMetaColor(object),label = F,
           split.by = "orig.ident",
           label.repel = T, pt.size = 0.5,label.size = 5, repel = T,no.legend = F,
           do.print = T,do.return = F,title = "Cell types in all 3 samples")

UMAPPlot.1(object, group.by = "orig.ident",cols = c("#ffa500","#0000ff","#008000"),label = T,
           label.repel = T, pt.size = 0.5,label.size = 5, repel = T,no.legend = T,
           do.print = T,do.return = F,title = "3 samples")

table(object$cell.types, object$orig.ident) %>% kable() %>%  kable_styling()
save(object, file = "data/Coloratoral_3_20191203.Rda")
