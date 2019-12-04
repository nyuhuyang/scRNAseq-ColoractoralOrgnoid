########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gplots)
library(MAST)
library(future)
source("../R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# change the current plan to access parallelization
plan("multiprocess", workers = 8)
plan()

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.numeric(slurm_arrayid)
print(paste0("slurm_arrayid=",args))

(load(file = paste0("data/Coloratoral_3_20191203.Rda")))
DefaultAssay(object)  = "SCT"
cell.types <- unique(sort(object$cell.types))
cell.type = cell.types[args]

Idents(object) = "cell.types"
sub_object <- subset(object, idents = cell.type)
Idents(sub_object) = "orig.ident"
system.time(markers <- FindAllMarkers.UMI(sub_object, 
                                          logfc.threshold = 0, only.pos = F,
                                          test.use = "MAST",
                                          min.cells.group = 1))
write.csv(markers,paste0(path,"KRPshS_",cell.type,"_markers_FC0.csv"))
