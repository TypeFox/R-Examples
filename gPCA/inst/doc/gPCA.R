### R code from vignette source 'gPCA.Rnw'

###################################################
### code chunk number 1: gPCA.Rnw:77-79 (eval = FALSE)
###################################################
## out<-gPCA.batchdetect(x=data,batch=batch,center=FALSE,
## 	scaleY=FALSE,filt=NULL,nperm=1000,seed=NULL)


###################################################
### code chunk number 2: gPCA.Rnw:90-94 (eval = FALSE)
###################################################
## gDist(out)
## PCplot(out,ug="guided",type="1v2")
## PCplot(out,ug="guided",type="comp",npcs=3)
## CumulativeVarPlot(out,ug="unguided",col="blue")


###################################################
### code chunk number 3: gPCA.Rnw:102-103
###################################################
library(gPCA)


###################################################
### code chunk number 4: gPCA.Rnw:105-112
###################################################
data(caseDat)
batch<-caseDat$batch
data<-caseDat$data

out<-gPCA.batchdetect(x=data,batch=batch,center=TRUE)
out$delta ; out$p.val
((out$varPCg1-out$varPCu1)/out$varPCg1)*100


###################################################
### code chunk number 5: gPCA.Rnw:120-121
###################################################
gDist(out)


###################################################
### code chunk number 6: gPCA.Rnw:132-134
###################################################
par(mai=c(0.8,0.8,0.1,0.1),cex=0.8)
PCplot(out,ug="guided",type="1v2")


###################################################
### code chunk number 7: gPCA.Rnw:143-145
###################################################
par(mai=c(0.65,0.65,0.1,0.1),cex=0.8)
PCplot(out,ug="guided",type="comp",npcs=3)


###################################################
### code chunk number 8: gPCA.Rnw:157-158
###################################################
sessionInfo()


