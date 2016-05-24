## ----setup1, include=FALSE-----------------------------------------------
require(knitr)
opts_chunk$set(fig.path='figure/beamer-',fig.align='center',fig.show='hold',size='footnotesize')

## ----setup2, include=FALSE------------------------------------------
options(width=70)  # make the printing fit on the page
set.seed(1121)     # make the results repeatable
stdt<-date()

## ----genes----------------------------------------------------------
genes <- c("BRCA1","BRCA2")

## ----snps-----------------------------------------------------------
snpInfo <- data.frame(chr=c(17,17,13,13),
                      pos=c(41211653, 41213996, 32890026,32890572),
                      rsid=c("rs8176273","rs8176265","rs9562605","rs1799943"),
                      stringsAsFactors=FALSE)

## ----loadpkg, size='tiny'-------------------------------------------
library(snplist)

## ----biomart1, eval=FALSE-------------------------------------------
#  geneInfo <- getBioMartData(genes)
#  geneInfo

## ----biomart2, echo=FALSE-------------------------------------------
geneInfo <- cbind(c('BRCA1','BRCA2'), c(17,13),c(41196312,32889611), 
		  c(41277500,32973805))
colnames(geneInfo) <- c('gene','chr','start','end')
geneInfo <- as.data.frame(geneInfo)
geneInfo

## ----genetbl--------------------------------------------------------
setGeneTable(geneInfo) 

## ----snptbl---------------------------------------------------------
setSNPTable(snpInfo)

## ----snpset1, size='tiny'-------------------------------------------
geneset <- makeGeneSet(margin=50000) 

## ----snpset2, tidy=FALSE, size='tiny'-------------------------------
chipSNPs <- c("rs0000000","rs8176273","rs9562605")
geneset2 <- makeGeneSet(annoInfo=chipSNPs, 
			annoTable='myChip') 

## ----plink----------------------------------------------------------
exportPLINKSet(geneset2, "mySNPset.set") 

## ----sesinf, echo=FALSE, results="asis"-----------------------------
toLatex(sessionInfo(), locale=FALSE)

## ----times, echo=FALSE----------------------------------------------
print(paste("Start Time",stdt))
print(paste("End Time  ",date()))

