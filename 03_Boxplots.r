## Boxplot 
library(ggplot2)
library(ggbeeswarm)

ggplot(tumor,aes(metastasis,meta.score))+geom_boxplot(outlier.shape=NA) + 
  scale_fill_brewer(palette="Dark2")+
  geom_point(aes(fill=factor(pretext_stage),size=1),shape = 21, alpha = .8, position = position_dodge2(width = .5))+
  theme_classic(base_size=18) +
  theme(legend.position = "right")+xlab("metastasis")+ggtitle("GSE131329-HB_tumors")


ggplot(tumor,aes(meta.cat,meta.score))+geom_boxplot(outlier.shape=NA) + 
  scale_fill_brewer(palette="Dark2")+
  geom_point(aes(fill=factor(pretext_stage),size=1),shape = 21, alpha = .8, position = position_dodge2(width = .5))+
  theme_classic(base_size=18) +
  theme(legend.position = "right")+ggtitle("GSE131329-HB_tumors")