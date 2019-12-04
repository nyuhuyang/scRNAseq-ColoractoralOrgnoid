library(SingleR)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/Coloratoral_data_4_20191008.Rda"))
mouse.rnaseq <-  MouseRNAseqData()
singler = SingleR(test = object_data, ref = mouse.rnaseq, 
                  labels = mouse.rnaseq$label.main)
save(singler,file="output/singler_T_Coloratoral_data_4_20191008")
  
