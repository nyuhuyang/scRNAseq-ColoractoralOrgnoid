library(Seurat)
library(magrittr)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggthemes)
source("../R/Seurat3_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
(load(file = "data/Coloratoral_4_20191008.Rda"))
(load(file="output/singler_T_Coloratoral_data_4_20191008.Rda"))

singler@nrows == ncol(object)
# if singler didn't find all cell labels
if(singler@nrows < ncol(object)) object = subset(object, cells = singler@rownames)
table(singler@rownames == colnames(object))

##############################
# add singleR label to Seurat
###############################

singlerDF = data.frame("labels" = singler$labels,
                       row.names = singler@rownames)
head(singlerDF)
table(singlerDF$labels) %>% kable() %>% kable_styling()
object <- AddMetaData(object = object,metadata = singlerDF)
object <- AddMetaColor(object = object, label= "labels", colors = Singler.colors)
Idents(object) <- "labels"

UMAPPlot.1(object, cols = ExtractMetaColor(object),label = T, label.repel = T,pt.size = 1,
         label.size = 4, repel = T,no.legend = T,do.print = T,
         title = "Major cell types")
save(object,file="data/Coloratoral_4_20191008.Rda")
##############################
# draw tsne plot
##############################
UMAPPlot.1(object, cols = ExtractMetaColor(object),split.by = "orig.ident",
           label = F, label.repel = F,pt.size = 1,border = T, ncol = 2,
           label.size = 4, repel = T,no.legend = T,do.print = T,
           title = "Major cell types",unique.name = "conditions")

