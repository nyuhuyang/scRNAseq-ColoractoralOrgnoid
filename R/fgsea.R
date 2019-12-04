########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(dplyr)
library(tidyr)
library(kableExtra)
library(magrittr)
library(gplots)
library(ggplot2)
library(fgsea)
library(tibble)
library(ggpubr)
library(ggsci)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 3.1.1 load data
# Load some results from Seurat
CSV_files <- list.files("Yang/DEGs")
cell.types = gsub("KRPshS_","",CSV_files)
(cell.types = gsub("_markers_FC0.csv","",cell.types))
hallmark <- gmtPathways("data/Intestinal_HALLMARK.gmt")
names(hallmark) = gsub("_"," ",names(hallmark))

hallmark %>% head() %>% lapply(head)

res <- list()
for(i in seq_along(cell.types)){
        res[[i]] = read.csv(file = paste0("Yang/DEGs/",CSV_files[i]),
                       row.names = 1, stringsAsFactors=F)
        res[[1]] = res[[1]][order(res[[1]]$avg_logFC, decreasing = T),]
        (clusters <- unique(res[[1]]$cluster))

        # Now, run the fgsea algorithm with 1000 permutations:
        fgseaRes = FgseaDotPlot(stats=res[[1]], pathways=hallmark, nperm=1000,padj = 0.25,pval = 0.05,
                     order.by = c("PRE","NES"),decreasing = F,
                     order = c("PRE","POST","RES"),
                     size = " -log10(pval)", fill = "NES",sample = cell.types[i], 
                     pathway.name = "Hallmark",rotate.x.text = F,
                     font.xtickslab=15, font.main=17, font.ytickslab = 10,
                     font.legend = list(size = 15),font.label = list(size = 15),
                     do.return = T,
                     save_path = paste0("Yang/GSEA/Dotplot_",cell.types[i],"_Hallmark_0_0.25_0.05.jpeg"),
                     width = 8,height = 8)
        write.csv(fgseaRes, file = paste0("Yang/GSEA/",cell.types[i],"_Hallmark.csv"))
}

