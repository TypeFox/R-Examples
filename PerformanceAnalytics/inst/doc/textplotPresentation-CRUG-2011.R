### R code from vignette source 'textplotPresentation-CRUG-2011.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Prelim
###################################################
library('PerformanceAnalytics')
data(managers)
#managers=read.csv("/home/peter/dev/R/managers.csv",row.names=1)
head(managers)
dim(managers)
colnames(managers)


###################################################
### code chunk number 2: Prelim
###################################################
manager.col = 1
peers.cols = c(2,3,4,5,6)
indexes.cols = c(7,8)
Rf.col = 10
peer.colorset=c("red", rep("darkorange", 2), rep("gray", 5))
ham1.downside = t(table.DownsideRisk(managers[,c(manager.col, 
indexes.cols, peers.cols)],Rf=.03/12))


###################################################
### code chunk number 3: ConstructTableEx
###################################################
ham1.downside


###################################################
### code chunk number 4: gplotstextplot
###################################################
library(gplots)
#args(gplots:::textplot)
gplots:::textplot(ham1.downside); box(col="lightblue")


###################################################
### code chunk number 5: Hmiscformat
###################################################
library(Hmisc)
args(format.df)
ham1.f.downside = format.df(ham1.downside, na.blank=TRUE, numeric.dollar = FALSE, cdec=rep(4,dim(ham1.downside)[2]))


###################################################
### code chunk number 6: Hmiscformat
###################################################
ham1.f.downside


###################################################
### code chunk number 7: PAtextplot
###################################################
require(PerformanceAnalytics)
args(PerformanceAnalytics:::textplot)


###################################################
### code chunk number 8: PAtextplot
###################################################
PerformanceAnalytics:::textplot(ham1.f.downside,  halign = "center", valign = "top", row.valign="center", col.rownames=peer.colorset,  mar = c(0,0,3,0)+0.1)
box(col="lightblue")


