library(readr)
library(pheatmap)
library(ggplot2)
library(hrbrthemes)
library(corrplot)

getwd()

#data from INCOR/BR
dnam <- read.csv("brs_dataset_incor_laval_norm.txt", sep = "\t")
head(dnam)
row.names(dnam)<- dnam[,1]
dnam<- dnam[,-1]
head(dnam)
dim(dnam)
dnam$Class <- factor(dnam$Class, levels=c("BrS", "Possible", "Control"))


mm<- data.matrix(dnam[,c(3:5)])
dim(mm)


cat_df = data.frame("Classification"=dnam$Class)
cat_df$Classification <- factor(cat_df$Classification, levels=c("BrS", "Possible", "Control"))
head(cat_df)
rownames(cat_df) = row.names(dnam)

pheatmap(log(t(mm)), cutree_cols = 3, cluster_rows = T, cluster_cols = F, annotation_col = cat_df, show_colnames =F,
         cellwidth = 15, cellheight = 18, fontsize = 8, 
         main= "Cluster autoantibodys log (O.D.) by Shanghay Classification Laval and INCOR samples", filename = "Heatmap_incor_laval.jpeg")


laval<- subset(dnam, id_2 == "Laval")
control_data<- subset(dnam, Class == "Control")
dnam_laval<- rbind(laval, control_data)
head(dnam_laval, 100)

mm<- data.matrix(dnam_laval[,c(3:5)])
dim(mm)

pheatmap(log(t(mm)), cutree_cols = 3, cluster_rows = T, cluster_cols = F, annotation_col = cat_df, show_colnames =F,
         cellwidth = 15, cellheight = 18, fontsize = 8, 
         main= "Cluster autoantibodys log (O.D.) by Shanghay Classification Laval samples", filename = "Heatmap_laval.jpeg")

incor<- subset(dnam, id_2 == "Incor")

mm<- data.matrix(incor[,c(3:5)])
dim(mm)

pheatmap(log(t(mm)), cutree_cols = 3, cluster_rows = T, cluster_cols = F, annotation_col = cat_df, show_colnames =F,
         cellwidth = 15, cellheight = 18, fontsize = 8, 
         main= "Cluster autoantibodys log (O.D.) by Shanghay Classification INCOR samples", filename = "Heatmap_INCOR.jpeg")