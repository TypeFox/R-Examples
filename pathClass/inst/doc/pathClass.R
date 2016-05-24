### R code from vignette source 'pathClass.rnw'

###################################################
### code chunk number 1: loadPackage
###################################################
library(pathClass)


###################################################
### code chunk number 2: loadALLpackage
###################################################
library(golubEsets)
data(Golub_Merge)


###################################################
### code chunk number 3: classLabels
###################################################
y <- pData(Golub_Merge)$ALL


###################################################
### code chunk number 4: loadDataMatrix
###################################################
x <- exprs(Golub_Merge)


###################################################
### code chunk number 5: dataMatrix
###################################################
x <- t(x)
dim(x)


###################################################
### code chunk number 6: hprd (eval = FALSE)
###################################################
## hprd <- read.hprd('BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt')


###################################################
### code chunk number 7: hprd2 (eval = FALSE)
###################################################
## data(adjacency.matrix)
## hprd <- adjacency.matrix


###################################################
### code chunk number 8: mapping
###################################################
ann <- annotation(Golub_Merge)
library(paste(ann,'db',sep='.'), character.only=TRUE)

graphIDs <- "REFSEQ"
rs <- get(paste(ann, graphIDs, sep=''))
refseq <- mget(featureNames(Golub_Merge), rs)
times <- sapply(refseq, length)
mapping <- data.frame(probesetID=rep(names(refseq), times=times),
                      graphID=unlist(refseq),
                      row.names=NULL,
                      stringsAsFactors=FALSE)
nas <- which(is.na(mapping[,'graphID']))
mapping <- mapping[-nas,]
mapping <- unique(mapping)
head(mapping)


###################################################
### code chunk number 9: matchMatrices (eval = FALSE)
###################################################
## matched <- matchMatrices(x=x, adjacency=hprd, mapping=mapping)


###################################################
### code chunk number 10: RRFE (eval = FALSE)
###################################################
## set.seed(12345)
## res.rrfe <- crossval(x,
##                      y,
##                      DEBUG=TRUE,
##                      theta.fit=fit.rrfe,
##                      folds=10,
##                      repeats=5,
##                      parallel=TRUE,
##                      Cs=10^(-3:3),
##                      mapping=mapping,
##                      Gsub=hprd,
##                      d=1/2)


###################################################
### code chunk number 11: RRFEallFeatures (eval = FALSE)
###################################################
## res.rrfe <- crossval(x,
##                      y,
##                      DEBUG=TRUE,
##                      theta.fit=fit.rrfe,
##                      folds=10,
##                      repeats=5,
##                      parallel=TRUE,
##                      Cs=10^(-3:3),
##                      useAllFeatures=TRUE,
##                      mapping=mapping,
##                      Gsub=hprd,
##                      d=1/2)


###################################################
### code chunk number 12: networkBasedSVM (eval = FALSE)
###################################################
## ad.list <- as.adjacencyList(matched$adjacency)
## 
## set.seed(12345)
## res.nBSVM <- crossval(matched$x,
##                       y,
##                       theta.fit=fit.networkBasedSVM,
##                       folds=10,
##                       repeats=5,
##                       DEBUG=TRUE,
##                       parallel=FALSE,
##                       adjacencyList=ad.list,
##                       lambdas=10^(-1:2),
##                       sd.cutoff=150)


###################################################
### code chunk number 13: graphSVM (eval = FALSE)
###################################################
## dk <- calc.diffusionKernel(L=matched$adjacency,
##                            is.adjacency=TRUE,
##                            beta=0)
## 
## set.seed(12345)
## res.gSVM <- crossval(matched$x,
##                      y,
##                      theta.fit=fit.graph.svm,
##                      folds=10,
##                      repeats=5,
##                      DEBUG=TRUE,
##                      parallel=FALSE,
##                      Cs=10^(-3:3),
##                      mapping=matched$mapping,
##                      diffusionKernel=dk)


###################################################
### code chunk number 14: plotting (eval = FALSE)
###################################################
## plot(res.rrfe, toFile=F)


###################################################
### code chunk number 15: benchmark (eval = FALSE)
###################################################
## cv.labels <- matrix(rep(y,5), ncol=5)
## pred.rrfe <- prediction(res.rrfe$cv, labels=cv.labels)
## auc.rrfe  <- round(mean(unlist(performance(pred.rrfe, 'auc')@y.values)),3)
## plot(performance(pred.rrfe, measure = "tpr", x.measure = "fpr"),
##      col='red',
##      main='Benchmark of the algorithms',
##      avg = "threshold")
## 
## pred.nBSVM <- prediction(res.nBSVM$cv, labels=cv.labels)
## auc.nBSVM  <- round(mean(unlist(performance(pred.nBSVM, 'auc')@y.values)),3)
## plot(performance(pred.nBSVM, measure = "tpr", x.measure = "fpr"),
##      add=TRUE,
##      col='blue',
##      avg = "threshold")
## 
## pred.gSVM <- prediction(res.gSVM$cv, labels=cv.labels)
## auc.gSVM  <- round(mean(unlist(performance(pred.gSVM, 'auc')@y.values)),3)
## plot(performance(pred.gSVM, measure = "tpr", x.measure = "fpr"),
##      add=TRUE,
##      col='green',
##      avg = "threshold")
## 
## legend('bottomright',
##        c(paste('RRFE (AUC=',auc.rrfe,')',sep=''),
##          paste('network based SVM (AUC=',auc.nBSVM,')',sep=''),
##          paste('graph SVM (AUC=',auc.gSVM,')',sep='')),
##        text.col=c('red','blue','green'),
##        col=c('red','blue','green'),
##        lty=1,
##        bty='n',
##        cex=1.3)
## 
## abline(b=1,a=0,col='gray')


###################################################
### code chunk number 16: extrFeatures (eval = FALSE)
###################################################
## extractFeatures(res.rrfe, toFile=T, fName='OurFeatures.csv')


