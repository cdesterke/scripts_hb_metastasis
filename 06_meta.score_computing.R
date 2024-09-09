### meta.score
data<-read.csv("input_data_meta.csv",h=T)

library(logitloop)
library(dplyr)


head(data)
data%>%filter(tissue=="tumor_tissue")->tumor

tumor%>%dplyr::rename(metastasis="distant_metastasis.ch1")->tumor

tumor$metastasis<-as.factor(tumor$metastasis)

tumor%>%select(2:42,metastasis)->small
df<-logitloop(small,outcome="metastasis")
df
plotcoef(df,nb=15,title="")
 df%>%filter(p.values<=0.1)%>%filter(beta>=0)->df
df<-df[1:2,]
id<-df$predictors
beta<-df$beta

equation<-paste(id,beta,sep="*",collapse=")+(")
equation<-paste("(",equation,")",sep="")
equation

tumor%>%mutate(meta.score=((DNMT3B*3.38905474467193)+(PFKFB4*2.3096192085631)))->tumor

library(cutpointr)


cp <- cutpointr(tumor, meta.score,  metastasis,method = maximize_metric, metric = sum_sens_spec)

plot(cp)
cp
tumor%>%mutate(meta.cat=ifelse(meta.score>=38.3052,"HIGH","low"))->tumor

chisq.test(table(tumor$meta.cat,tumor$clinical_course))

tumor$meta.cat<-as.factor(tumor$meta.cat)
tumor$meta.cat<-relevel(tumor$meta.cat,ref="low")