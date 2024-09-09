## ROC curve on positive coefficients
library(pROC)
library(Epi)
data<-as.data.frame(cbind(X,y))
ROC(form = meta~SOAT2+DDAH2+PYCR1+SOD3+PFKFB4+DNMT3B+ISYNA1+FKBP10+PKM+NT5DC2+CHST10+
ENO2+GPX7+GSTP1+PAPSS1+HK2+PFKM+P4HA2 , plot="ROC", data=data)