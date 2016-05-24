### R code from vignette source 'amap.Rnw'

###################################################
### code chunk number 1: amap.Rnw:32-36
###################################################
library(amap)
data(USArrests)
h = hcluster(USArrests)
plot(h)


###################################################
### code chunk number 2: amap.Rnw:39-42
###################################################
heatmap(as.matrix(USArrests),
        hclustfun=hcluster,
        distfun=function(u){u})


###################################################
### code chunk number 3: amap.Rnw:45-46
###################################################
h = hcluster(USArrests,nbproc=4)   


###################################################
### code chunk number 4: amap.Rnw:49-50
###################################################
Kmeans(USArrests,centers=3,method="correlation")


###################################################
### code chunk number 5: amap.Rnw:56-59
###################################################
data(lubisch)
lubisch <- lubisch[,-c(1,8)]
varrob(scale(lubisch),h=1)


###################################################
### code chunk number 6: amap.Rnw:62-64
###################################################
p <- acpgen(lubisch,h1=1,h2=1/sqrt(2))
plot(p)


###################################################
### code chunk number 7: amap.Rnw:67-69
###################################################
p <- acprob(lubisch,h=4)
plot(p)


