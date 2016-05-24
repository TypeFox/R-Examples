### R code from vignette source 'metaRNASeq.Rnw'

###################################################
### code chunk number 1: metaRNASeq.Rnw:42-43
###################################################
options(width=60)


###################################################
### code chunk number 2: loadparameters
###################################################
library(metaRNASeq)
library(DESeq2)
data(param)
dim(param)
data(dispFuncs)


###################################################
### code chunk number 3: simulateData
###################################################
set.seed(123)
matsim <- sim.function(param = param, dispFuncs = dispFuncs)
sim.conds <- colnames(matsim)
rownames(matsim) <- paste("tag", 1:dim(matsim)[1],sep="")
dim(matsim)


###################################################
### code chunk number 4: extractindivstudy
###################################################
colnames(matsim)
simstudy1 <- extractfromsim(matsim,"study1")
head(simstudy1$study)
simstudy1$pheno
simstudy2 <- extractfromsim(matsim,"study2")


###################################################
### code chunk number 5: DESeq2.indivanalysis
###################################################
dds1 <- DESeqDataSetFromMatrix(countData = simstudy1$study,
 colData = simstudy1$pheno,design = ~ condition)
res1 <- results(DESeq(dds1))
dds2 <- DESeqDataSetFromMatrix(countData = simstudy2$study, 
  colData = simstudy2$pheno,design = ~ condition)
res2 <- results(DESeq(dds2))


###################################################
### code chunk number 6: storepvalandFC
###################################################
rawpval <- list("pval1"=res1[["pvalue"]],"pval2"=res2[["pvalue"]])
FC <- list("FC1"=res1[["log2FoldChange"]],"FC2"=res2[["log2FoldChange"]]) 


###################################################
### code chunk number 7: storeadjpval
###################################################
adjpval <- list("adjpval1"=res1[["padj"]],"adjpval2"=res2[["padj"]])
DE <- mapply(adjpval, FUN=function(x) ifelse(x <= 0.05, 1, 0))
colnames(DE)=c("DEstudy1","DEstudy2")


###################################################
### code chunk number 8: pvalDESeq2hist
###################################################
par(mfrow = c(1,2))
hist(rawpval[[1]], breaks=100, col="grey", main="Study 1", xlab="Raw p-values")
hist(rawpval[[2]], breaks=100, col="grey", main="Study 2", xlab="Raw p-values")


###################################################
### code chunk number 9: DESeq1study
###################################################
library(DESeq)
library(HTSFilter)
resDESeq1study <- function(studyname, alldata, cond1totest="cond1",
    cond2totest="cond2", fitType = "parametric") {
  study <- alldata[,grep(studyname,colnames(alldata))]
  studyconds <- gsub(studyname,"",colnames(study))
  colnames(study) <- paste(studyconds,1:dim(study)[2],sep=".")
  cds <- newCountDataSet(study, studyconds)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds, method="pooled", fitType=fitType)
  ## Filter using Jaccard index for each study
  filter <- HTSFilter(cds, plot=FALSE)
  cds.filter <- filter$filteredData
  on.index <- which(filter$on == 1)
  cat("# genes passing filter", studyname, ":", dim(cds.filter)[1], "\n")
  res <- as.data.frame(matrix(NA, nrow = nrow(cds), ncol=ncol(cds)))
  nbT <- nbinomTest(cds.filter, cond1totest, cond2totest)
  colnames(res) <- colnames(nbT)
  res[on.index,] <- nbT
  res
}


###################################################
### code chunk number 10: DESeq2studies
###################################################
studies <- c("study1", "study2")
resDESeq <- lapply(studies, 
  FUN=function(x) resDESeq1study(x, alldata=matsim))


###################################################
### code chunk number 11: pvalDE
###################################################
rawpval <- lapply(resDESeq, FUN=function(res) res$pval)
adjpval <- lapply(resDESeq, FUN=function(res) res$padj)
DE <- mapply(adjpval, FUN=function(x) ifelse(x <= 0.05, 1, 0))
colnames(DE)=paste("DE",studies,sep=".")


###################################################
### code chunk number 12: pvalDEhist
###################################################
par(mfrow = c(1,2))
hist(rawpval[[1]], breaks=100, col="grey", main="Study 1", 
  xlab="Raw p-values")
hist(rawpval[[2]], breaks=100, col="grey", main="Study 2", 
  xlab="Raw p-values")


###################################################
### code chunk number 13: pvalfishcomb
###################################################
fishcomb <- fishercomb(rawpval, BHth = 0.05)
hist(fishcomb$rawpval, breaks=100, col="grey", main="Fisher method",
  xlab = "Raw p-values (meta-analysis)")


###################################################
### code chunk number 14: pvalinvnorm
###################################################
invnormcomb <- invnorm(rawpval,nrep=c(8,8), BHth = 0.05)   
hist(invnormcomb$rawpval, breaks=100, col="grey", 
  main="Inverse normal method",
  xlab = "Raw p-values (meta-analysis)")    


###################################################
### code chunk number 15: tabDE
###################################################
DEresults <- data.frame(DE, 
  "DE.fishercomb"=ifelse(fishcomb$adjpval<=0.05,1,0),
  "DE.invnorm"=ifelse(invnormcomb$adjpval<=0.05,1,0))
head(DEresults)


###################################################
### code chunk number 16: checkDESeq2
###################################################
signsFC <- mapply(FC, FUN=function(x) sign(x))
sumsigns <- apply(signsFC,1,sum) 
commonsgnFC <- ifelse(abs(sumsigns)==dim(signsFC)[2], sign(sumsigns),0)  


###################################################
### code chunk number 17: filterconflicts
###################################################
unionDE <- unique(c(fishcomb$DEindices,invnormcomb$DEindices))
FC.selecDE <- data.frame(DEresults[unionDE,],do.call(cbind,FC)[unionDE,],
  signFC=commonsgnFC[unionDE], DE=param$DE[unionDE])
keepDE <- FC.selecDE[which(abs(FC.selecDE$signFC)==1),]
conflictDE <- FC.selecDE[which(FC.selecDE$signFC == 0),]
dim(FC.selecDE)
dim(keepDE)
dim(conflictDE)
head(keepDE)


###################################################
### code chunk number 18: filtercheckcache
###################################################
nbtrueconflicts=as.vector(table(conflictDE$DE)[2])


###################################################
### code chunk number 19: filtercheck
###################################################
table(conflictDE$DE)


###################################################
### code chunk number 20: calcul
###################################################
fishcomb_de <- rownames(keepDE)[which(keepDE[,"DE.fishercomb"]==1)] 
invnorm_de <- rownames(keepDE)[which(keepDE[,"DE.invnorm"]==1)] 
indstudy_de <- list(rownames(keepDE)[which(keepDE[,"DE.study1"]==1)], 
                    rownames(keepDE)[which(keepDE[,"DE.study2"]==1)])

IDD.IRR(fishcomb_de,indstudy_de)
IDD.IRR(invnorm_de ,indstudy_de)


###################################################
### code chunk number 21: calcul2
###################################################
x=IDD.IRR(fishcomb_de,indstudy_de)
y=IDD.IRR(invnorm_de ,indstudy_de)


###################################################
### code chunk number 22: venndiagram
###################################################
library(VennDiagram)
venn.plot<-venn.diagram(x = list(study1=which(keepDE[,"DE.study1"]==1),
                                 study2=which(keepDE[,"DE.study2"]==1),
                                 fisher=which(keepDE[,"DE.fishercomb"]==1),
                                 invnorm=which(keepDE[,"DE.invnorm"]==1)),
                        filename = NULL, col = "black",
                        fill = c("blue", "red", "purple","green"),
                        margin=0.05, alpha = 0.6)
jpeg("venn_jpeg.jpg");
grid.draw(venn.plot);
dev.off();


###################################################
### code chunk number 23: sessionInfo
###################################################
sessionInfo()


