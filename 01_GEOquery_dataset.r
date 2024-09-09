## get public dataset GSE131329
library(GEOquery)
# download data GEO GSE131329
gse <- getGEO("GSE131329", GSEMatrix = TRUE)

# if dataset has several data you can access like
set <- gse[[1]]

# extract expression
data<-exprs(set)

# extract annotation
features<-fData(set)

# extract phenotype
annot <- pData(set)


# export data
all<-merge(features,data,by="row.names")

write.csv(all,file="data.csv",row.names=F)
write.csv(annot,file="pheno.csv",row.names=T)