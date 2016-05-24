### R code from vignette source 'clValid.Rnw'

###################################################
### code chunk number 1: prelim
###################################################
options(prompt="R> ")
library("clValid")


###################################################
### code chunk number 2: mouse
###################################################
data(mouse)


###################################################
### code chunk number 3: internal
###################################################
express <- mouse[,c("M1","M2","M3","NC1","NC2","NC3")]
rownames(express) <- mouse$ID
intern <- clValid(express, 2:6, clMethods=c("hierarchical","kmeans","pam"),
                  validation="internal")


###################################################
### code chunk number 4: clValid.Rnw:659-660
###################################################
summary(intern)


###################################################
### code chunk number 5: internPlot
###################################################
op <- par(no.readonly=TRUE)
par(mfrow=c(2,2),mar=c(4,4,3,1))
plot(intern, legend=FALSE)
plot(nClusters(intern),measures(intern,"Dunn")[,,1],type="n",axes=F,
        xlab="",ylab="")
legend("center", clusterMethods(intern), col=1:9, lty=1:9, pch=paste(1:9))
par(op)


###################################################
### code chunk number 6: clValid.Rnw:686-687
###################################################
op <- par(no.readonly=TRUE)
par(mfrow=c(2,2),mar=c(4,4,3,1))
plot(intern, legend=FALSE)
plot(nClusters(intern),measures(intern,"Dunn")[,,1],type="n",axes=F,
        xlab="",ylab="")
legend("center", clusterMethods(intern), col=1:9, lty=1:9, pch=paste(1:9))
par(op)


###################################################
### code chunk number 7: stability
###################################################
stab <- clValid(express, 2:6, clMethods=c("hierarchical","kmeans","pam"),
                validation="stability")


###################################################
### code chunk number 8: optimal
###################################################
optimalScores(stab)


###################################################
### code chunk number 9: stabPlot
###################################################
par(mfrow=c(2,2),mar=c(4,4,3,1))
plot(stab, measure=c("APN","AD","ADM"),legend=FALSE)
plot(nClusters(stab),measures(stab,"APN")[,,1],type="n",axes=F,
        xlab="",ylab="")
legend("center", clusterMethods(stab), col=1:9, lty=1:9, pch=paste(1:9))
par(op)


###################################################
### code chunk number 10: clValid.Rnw:749-750
###################################################
par(mfrow=c(2,2),mar=c(4,4,3,1))
plot(stab, measure=c("APN","AD","ADM"),legend=FALSE)
plot(nClusters(stab),measures(stab,"APN")[,,1],type="n",axes=F,
        xlab="",ylab="")
legend("center", clusterMethods(stab), col=1:9, lty=1:9, pch=paste(1:9))
par(op)


###################################################
### code chunk number 11: biological
###################################################
fc <- tapply(rownames(express),mouse$FC, c)
fc <- fc[!names(fc)%in%c("EST","Unknown")]
bio <- clValid(express, 2:6, clMethods=c("hierarchical","kmeans","pam"),
               validation="biological", annotation=fc)


###################################################
### code chunk number 12: readExternalBiological
###################################################
fc2 <- readAnnotationFile("fc.csv")
## bio.fc2 <- clValid(express, 2:6, clMethods=c("hierarchical","kmeans","pam"),
##                    validation="biological", annotation=fc2)
## all.equal(measures(bio), measures(bio.fc2))


###################################################
### code chunk number 13: bioOptimal
###################################################
optimalScores(bio)


###################################################
### code chunk number 14: BHI
###################################################
plot(bio, measure="BHI", legendLoc="topleft")


###################################################
### code chunk number 15: BSI
###################################################
plot(bio, measure="BSI")


###################################################
### code chunk number 16: clValid.Rnw:846-847
###################################################
plot(bio, measure="BHI", legendLoc="topleft")


###################################################
### code chunk number 17: clValid.Rnw:856-857
###################################################
plot(bio, measure="BSI")


###################################################
### code chunk number 18: bioc
###################################################
if(require("Biobase") && require("annotate") && require("GO.db") && require("moe430a.db")) {
  ## Need to know which affy chip was used in experiment
  ## affymetrix murine genome 430a genechip arrays
  bio2 <- clValid(express, 2:6, clMethods=c("hierarchical","kmeans","pam"),
                  validation="biological",annotation="moe430a.db",GOcategory="all")
}



###################################################
### code chunk number 19: clValid.Rnw:904-905
###################################################
if(exists("bio2")) optimalScores(bio2)


###################################################
### code chunk number 20: BHI2
###################################################
if(exists("bio2")) plot(bio2, measure="BHI", legendLoc="topleft")


###################################################
### code chunk number 21: BSI2
###################################################
if(exists("bio2")) plot(bio2, measure="BSI")


###################################################
### code chunk number 22: clValid.Rnw:926-927
###################################################
if(exists("bio2")) plot(bio2, measure="BHI", legendLoc="topleft")


###################################################
### code chunk number 23: clValid.Rnw:936-937
###################################################
if(exists("bio2")) plot(bio2, measure="BSI")


###################################################
### code chunk number 24: biocDE (eval = FALSE)
###################################################
## if(require("Biobase") && require("annotate") && require("GO.db") && require("moe430a.db")) {
##   ## Need to know which affy chip was used in experiment
##   ## affymetrix murine genome 430a genechip arrays
##   bio2DE <- clValid(express, 2:6, clMethods=c("hierarchical","kmeans","pam"),
##                   validation="biological",annotation="moe430a.db",GOcategory="all",
##                   dropEvidence="IEA")
##   optimalScores(bio2DE)
## }


###################################################
### code chunk number 25: RankAggreg
###################################################
result <- clValid(express, 4:6, clMethods=c("hierarchical","kmeans","pam"), 
                  validation=c("internal","stability"))
res <- getRanksWeights(result)


###################################################
### code chunk number 26: ranks
###################################################
print(res$ranks[,1:3], quote=FALSE)


###################################################
### code chunk number 27: RankAggreg
###################################################
if(require("RankAggreg")) {
  CEWS <- RankAggreg(x=res$ranks, k=5, weights=res$weights, seed=123, verbose=FALSE)
  CEWS
}


###################################################
### code chunk number 28: RankAggFig
###################################################
plot(CEWS)


###################################################
### code chunk number 29: RAFig
###################################################
plot(CEWS)


###################################################
### code chunk number 30: hierarchical
###################################################
hc <- clusters(bio,"hierarchical")


###################################################
### code chunk number 31: hplot
###################################################
mfc <- factor(mouse$FC, labels=c("Re","EST","GD","KP","Met","Mis","St","TF","U"))
tf.gd <- ifelse(mfc%in%c("GD","TF"),levels(mfc)[mfc],"")
plot(hc, labels=tf.gd, cex=0.7, hang=-1, main="Mouse Cluster Dendrogram")


###################################################
### code chunk number 32: clValid.Rnw:1093-1094
###################################################
mfc <- factor(mouse$FC, labels=c("Re","EST","GD","KP","Met","Mis","St","TF","U"))
tf.gd <- ifelse(mfc%in%c("GD","TF"),levels(mfc)[mfc],"")
plot(hc, labels=tf.gd, cex=0.7, hang=-1, main="Mouse Cluster Dendrogram")


###################################################
### code chunk number 33: twoClusters
###################################################
two <- cutree(hc,2)
xtabs(~mouse$FC + two)


