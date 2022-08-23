library(readr)

###HEATMAPS

library(pheatmap)
getwd()

#####data from INCOR/BR + Laval normalized by |O.D. value - means control|
dnam <- read.csv("brs_dataset_incor_laval_norm.txt", sep = "\t")
head(dnam)
row.names(dnam)<- dnam[,1]
dnam<- dnam[,-1]
head(dnam)
dim(dnam)
dnam$Class <- factor(dnam$Class, levels=c("BrS", "Possible", "Control"))

#####transfor in data matrix
mm<- data.matrix(dnam[,c(3:5)])
dim(mm)

#####creat vector classification
cat_df = data.frame("Classification"=dnam$Class)
cat_df$Classification <- factor(cat_df$Classification, levels=c("BrS", "Possible", "Control"))
head(cat_df)
rownames(cat_df) = row.names(dnam)

#####heatmap for all dataset
pheatmap(log(t(mm)), cutree_cols = 3, cluster_rows = T, cluster_cols = F, annotation_col = cat_df, show_colnames =F,
         cellwidth = 15, cellheight = 18, fontsize = 8, 
         main= "Cluster autoantibodys log (O.D.) by Shanghay Classification Laval and INCOR samples", filename = "Heatmap_incor_laval.jpeg")

#####subset data only Laval
laval<- subset(dnam, id_2 == "Laval")
control_data<- subset(dnam, Class == "Control")
dnam_laval<- rbind(laval, control_data)
head(dnam_laval)
mm<- data.matrix(dnam_laval[,c(3:5)])
dim(mm)

pheatmap(log(t(mm)), cutree_cols = 3, cluster_rows = T, cluster_cols = F, annotation_col = cat_df, show_colnames =F,
         cellwidth = 15, cellheight = 18, fontsize = 8, 
         main= "Cluster autoantibodys log (O.D.) by Shanghay Classification Laval samples", filename = "Heatmap_laval.jpeg")

#####subset data only INCOR
incor<- subset(dnam, id_2 == "Incor")
mm<- data.matrix(incor[,c(3:5)])
dim(mm)

pheatmap(log(t(mm)), cutree_cols = 3, cluster_rows = T, cluster_cols = F, annotation_col = cat_df, show_colnames =F,
         cellwidth = 15, cellheight = 18, fontsize = 8, 
         main= "Cluster autoantibodys log (O.D.) by Shanghay Classification INCOR samples", filename = "Heatmap_INCOR.jpeg")

################################################################################

### Correlation plot

library(corrplot)

#####data from INCOR/BR + Laval normalized by |O.D. value - means control|
dnam <- read.csv("brs_dataset_incor_laval_norm.txt", sep = "\t")
row.names(dnam)<- dnam[,1]
dnam<- dnam[,-1]
dnam$Class <- factor(dnam$Class, levels=c("BrS", "Possible", "Control"))
head(dnam)
dim(dnam)
#####transfor in data matrix
mm<- data.matrix(dnam[,c(3:6)])
dim(mm)
M<-cor(mm)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
png(height=800, width=800, file="Correlation O.D. All Dataset.png", type = "cairo")
corrplot(M, method="circle", col=col(100),  
         diag=FALSE,
         type="lower", order="FPC", 
         addCoef.col = "black",
         mar=c(0,0,1,0),
         tl.srt = 45,
         number.cex= 2.5,
         tl.cex = 2.3,
         cl.cex = 2, number.digits = 2)
dev.off()

#####subset data only Laval
laval<- subset(dnam, id_2 == "Laval")
control_data<- subset(dnam, Class == "Control")
dnam_laval<- rbind(laval, control_data)
head(dnam_laval)
#####transfor in data matrix
mm<- data.matrix(dnam_laval[,c(3:6)])
dim(mm)
M<-cor(mm)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
png(height=800, width=800, file="Correlation O.D. Laval Dataset.png", type = "cairo")
corrplot(M, method="circle", col=col(100),  
         diag=FALSE,
         type="lower", order="FPC", 
         addCoef.col = "black",
         mar=c(0,0,1,0),
         tl.srt = 45,
         number.cex= 2.5,
         tl.cex = 2.3,
         cl.cex = 2, number.digits = 3)
dev.off()

#####subset data only Incor
incor<- subset(dnam, id_2 == "Incor")
#####transfor in data matrix
mm<- data.matrix(incor[,c(3:6)])
dim(mm)
M<-cor(mm)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
png(height=800, width=800, file="Correlation O.D. INCOR Dataset.png", type = "cairo")
corrplot(M, method="circle", col=col(100),  
         diag=FALSE,
         type="lower", order="FPC", 
         addCoef.col = "black",
         mar=c(0,0,1,0),
         tl.srt = 45,
         number.cex= 2.5,
         tl.cex = 2.3,
         cl.cex = 2, number.digits = 3)
dev.off()

################################################################################


###PCA

pca = prcomp(t(mm))
class(pca)
str(pca)

pca$x
plot(pca$x)
text(pca$x, colnames(mm), col ='red')








install.packages("cutpointr")
install.packages("pROC")
library(pROC)
library(cutpointr)


dnam <- read.csv("brs_dataset_o.txt", sep = "\t")
row.names(dnam)<- dnam[,1]
dnam<- dnam[,-1]
head(dnam)


rocobj <- roc(dnam$response, dnam$actPosNeg)

rocobj
ggroc(rocobj)


abs_d_ppv_npv(19, 1, 2, 1)
accuracy(19, 1, 2, 1)




cp <- cutpointr(dnam$response , dnam$Actin , dnam$response, pos_class = "yes", neg_class = "no", direction = ">=",
                method = maximize_metric, metric = sum_sens_spec, na.rm=T)
summary(cp)



cat_df


#ROC Curve

library(plotROC)

dnam2 <- read.csv("brs_dataset_roc.txt", sep = "\t")
head(dnam2)

basicplot<- ggplot(dnam2, aes(d = D, m = M, color = Biomarker)) + geom_roc(n.cuts = 0) 
basicplot + style_roc()
  
annotate("text", x = .75, y = .25, 
         label = paste("AUC =", round(calc_auc(basicplot)$AUC, 4)))


