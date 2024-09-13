## C1/C2 signature classifier

data<-read.csv("matrix.csv", h=T)

data<-data[complete.cases(data$gene),]

library(transpipe15)
data<-filtermatrix(data)

pheno<-read.csv("input_data_meta.csv",h=T,row.names=1)

all(colnames(data)==row.names(pheno))

library(dplyr)
pheno%>%dplyr::rename(metastasis="distant_metastasis.ch1")->pheno


cairo<-c("GHR","APCS","C1S","AQP9","CYP2E1","APOC4","HDP","NLE","RPL10A","E2F5","BUB1","DLG7","IGSF1","AFP","DUSP9","ALDH2")
length(cairo)

inter<-intersect(cairo,row.names(data))
length(inter)

small<-data[row.names(data)%in%inter,]

pheno%>%filter(tissue!="noncancerous_liver_tissue")->pheno_tumor

cairo_tumor<-small[,row.names(pheno_tumor)]

trans<-as.data.frame(t(cairo_tumor))



##kmeans
tot.withinss <- vector(mode="character", length=10)
for (i in 1:10){
  Cluster <- kmeans(trans, center=i, nstart=20)
  tot.withinss[i] <- Cluster$tot.withinss
}
Let’s visualize it.

plot(1:10, tot.withinss, type="b", pch=19)

set.seed(101)
Cluster <- kmeans(trans, center=2, nstart=20)
Cluster

prop.table(table(pheno_tumor$meta.cat,Cluster$cluster))*100


tumor<-read.csv("tumor.csv",h=T)

all(row.names(pheno_tumor)== (tumor$id))

pheno_tumor$meta.cat<-tumor$meta.cat
pheno_tumor$meta.score<-tumor$meta.score


df<-table(pheno_tumor$meta.cat,Cluster$cluster)

pheno_tumor$kmeans<-Cluster$cluster
pheno_tumor$kmeans<-as.factor(pheno_tumor$kmeans)
pcatrans(cairo_tumor,pheno_tumor,group="cairo",pal="Dark2",names=F,alpha=1)




pheno_tumor%>%select(cairo,meta.cat,metastasis)->annot

bestheat(cairo_tumor,annot,scale="none",font=12)

pheno

pheno_tumor%>%mutate(cairo=case_when(kmeans=="1"~"C2",kmeans=="2"~"C1"))->pheno_tumor



library(cluster)
clusplot(trans, Cluster$cluster, color=T, shade=F, labels=0, lines=0)






library(caret)
cm <- confusionMatrix(data = pheno$clusters, reference = pheno$consensus)

pheno$consensus<-as.factor(pheno$consensus)

table(pheno$clusters,pheno$consensus)

draw_confusion_matrix(cm)





draw_confusion_matrix <- function(cm) {

  total <- sum(cm$table)
  res <- as.numeric(cm$table)

  # Generate color gradients. Palettes come from RColorBrewer.
  greenPalette <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
  redPalette <- c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
  getColor <- function (greenOrRed = "green", amount = 0) {
    if (amount == 0)
      return("#FFFFFF")
    palette <- greenPalette
    if (greenOrRed == "red")
      palette <- redPalette
    colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
  }

  # set the basic layout
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)

  # create the matrix 
  classes = colnames(cm$table)
  rect(150, 430, 240, 370, col=getColor("green", res[1]))
  text(195, 435, classes[1], cex=1.2)
  rect(250, 430, 340, 370, col=getColor("red", res[3]))
  text(295, 435, classes[2], cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=getColor("red", res[2]))
  rect(250, 305, 340, 365, col=getColor("green", res[4]))
  text(140, 400, classes[1], cex=1.2, srt=90)
  text(140, 335, classes[2], cex=1.2, srt=90)

  # add in the cm results
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')

  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)

  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}


library(caret)


pheno_tumor%>%mutate(cairo_code=case_when(cairo=="C2"~"1",cairo=="C1"~"0"))->pheno_tumor
pheno_tumor%>%mutate(metastasis_code=case_when(metastasis=="None"~"0",metastasis=="positive"~"1"))->pheno_tumor
pheno_tumor$metastasis_code<-as.factor(pheno_tumor$metastasis_code)
pheno_tumor$cairo_code<-as.factor(pheno_tumor$cairo_code)
cm <- confusionMatrix(data = pheno_tumor$cairo_code, reference = pheno_tumor$metastasis_code)


draw_confusion_matrix(cm)

pheno_tumor%>%mutate(metacat_code=case_when(meta.cat=="low"~"0",meta.cat=="HIGH"~"1"))->pheno_tumor
pheno_tumor$metacat_code<-as.factor(pheno_tumor$metacat_code)
cm <- confusionMatrix(data = pheno_tumor$metacat_code, reference = pheno_tumor$metastasis_code)


draw_confusion_matrix(cm)


write.table(pheno_tumor,file="pheno_tumor.tsv",row.names=T,sep="\t")

save(cairo_tumor,file="cairo.rda")

















