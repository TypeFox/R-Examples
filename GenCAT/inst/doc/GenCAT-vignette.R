## ---- message=FALSE------------------------------------------------------
library(GenCAT)

## ------------------------------------------------------------------------
data('CardioData')
data('coords')

print(head(coords))
print(head(CardioData))

CardioMapped<-map2class(coords, CardioData, extend.boundary = 5000)

print(head(CardioMapped))

## ------------------------------------------------------------------------
data('geno')

genoData<-geno$genotypes
snpInfo<-geno$map
colnames(snpInfo)<-c('chr', 'SNP', 'gen.dist', 'position', 'A1', 'A2')

GenCATtest <- GenCAT(CardioMapped, genoData=genoData, snpInfo = snpInfo)

## ------------------------------------------------------------------------
print(str(GenCATtest))

## ---- fig.align='center', fig.width = 10, fig.height=6-------------------
GenCAT_manhattan(GenCATtest, sigThresh = (0.05/nrow(GenCATtest$GenCAT)), 
highlightPosi = TRUE, labelPosi = TRUE)

