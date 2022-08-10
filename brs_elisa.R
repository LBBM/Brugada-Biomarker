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







ggplot(dnam, aes(x=actin, y=kera)) + 
  geom_point( color="#69b3a2") +
  theme_ipsum()
warnings()

ggplot(dnam, aes(x=log(actin), y=log(kera))) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  geom_text(x = -1, y = 0.7, label = eq(dnam$actin,dnam$kera), parse = TRUE)+
  theme_ipsum()


eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}




#data from LAVAL/CA

dnam <- read.csv("brs_dataset_laval.txt", sep = "\t")

row.names(dnam)<- dnam[,1]
dnam<- dnam[,-1]
head(dnam)
dim(dnam)
dnam$Class <- factor(dnam$Class, levels=c("Control", "Non-BrS", "Possible", "BrS"))
summary(dnam)

mm<- data.matrix(dnam[,c(1:4,12)])
dim(mm)
mm
cor(mm)

library(pheatmap)

cat_df = data.frame("Mutation" = dnam$mut, "Shangai"=dnam$Class)
rownames(cat_df) = row.names(dnam)
#all samples
pheatmap(log(t(mm+1)), cutree_cols = 5, cluster_rows = F, annotation_col = cat_df)






#only


pheatmap::pheatmap(mm, border_color = NA)

?pheatmap

cor(mm)

pheatmap(cor(mm), border_color = NA)


cor(mm[,1],mm[,2])
dim(mm)
nrow(mm)
ncol(mm)
plot(mm)

plot(mm[,1],mm[,3])


heatmap(mm)


heatmap(log(mm))


#PCA

pca = prcomp(t(mm))
class(pca)
str(pca)

pca$x
plot(pca$x)
text(pca$x, colnames(mm), col ='red')


#BloxPlot


library(ggstatsplot)

#dnam_brs<- subset(dnam, Class == "BrS" | Class == "Control")
#dnam_brs

dnam$classif <- factor(dnam$classif, levels=c("Control", "Possible", "BrS"))

plt <- ggbetweenstats(
  data = dnam,
  x = classif,
  y =conn, type = "np"
)
plt

plt <- plt + 
  # Add labels and title
  labs(
    x = "Shangai Score",
    y = "O.D. (Connexin)",
    #title = "Distribution of bill length across penguins species"
  ) + 
  # Customizations
  theme(
    # This is the new default font in the plot
    text = element_text(family = "Roboto", size = 8, color = "black"),
    plot.title = element_text(
      family = "Lobster Two", 
      size = 20,
      face = "bold",
      color = "#2a475e"
    ),
    # Statistical annotations below the main title
    plot.subtitle = element_text(
      family = "Roboto", 
      size = 15, 
      face = "bold",
      color="#1b2838"
    ),
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12)
  )

plt <- plt  +
  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
    plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
  )

plt


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


#Matrix corr plot
dnam <- read.csv("brs_dataset_corr.txt", sep = "\t")
head(dnam)
dim(dnam)
mm<- data.matrix(dnam)
dim(mm)

M<-cor(mm)
corrplot(M, method="circle")
