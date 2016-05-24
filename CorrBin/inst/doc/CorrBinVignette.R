### R code from vignette source 'CorrBinVignette.Rnw'

###################################################
### code chunk number 1: Misc
###################################################
options(useFancyQuotes=FALSE) 


###################################################
### code chunk number 2: LibLoad
###################################################
  library(CorrBin)
  library(lattice)


###################################################
### code chunk number 3: Intro
###################################################
 options(width=100)
 ps.options(colormodel="rgb")


###################################################
### code chunk number 4: DataLoad
###################################################
 sh <- read.CBData("ShellTox.txt", with.freq=TRUE)
 levels(sh$Trt) <- c("Control","Low","Medium", "High")
 str(sh)


###################################################
### code chunk number 5: MCtest
###################################################
  mc.test.chisq(sh)


###################################################
### code chunk number 6: MCest
###################################################
  sh.mc <- mc.est(sh)


###################################################
### code chunk number 7: MCestfig
###################################################
  print(xyplot(Prob~NResp|factor(ClusterSize), groups=Trt, data=sh.mc, subset=ClusterSize>0 & ClusterSize<13, 
    type="l", as.table=TRUE, auto.key=list(columns=4, lines=TRUE, points=FALSE),
    xlab="Number of responses", ylab="P(R=r|N=n)"))


###################################################
### code chunk number 8: MCestfig2 (eval = FALSE)
###################################################
##   print(xyplot(Prob~NResp|factor(ClusterSize), groups=Trt, data=sh.mc, subset=ClusterSize>0 & ClusterSize<13, 
##     type="l", as.table=TRUE, auto.key=list(columns=4, lines=TRUE, points=FALSE),
##     xlab="Number of responses", ylab="P(R=r|N=n)"))


###################################################
### code chunk number 9: MCsurvfig
###################################################
  panel.cumsum <- function(x,y,...){
    x.ord <- order(x)
    panel.xyplot(x[x.ord], cumsum(y[x.ord]), ...)}

   print(xyplot(Prob~NResp|factor(ClusterSize), groups=Trt, data=sh.mc, 
       subset=ClusterSize>0&ClusterSize<13, type="s",
       panel=panel.superpose, panel.groups=panel.cumsum,
       as.table=T, auto.key=list(columns=4, lines=T, points=F),
			 xlab="Number of responses", ylab="Cumulative Probability R(R>=r|N=n)",
			 ylim=c(0,1.1)))


###################################################
### code chunk number 10: MCsurvfig2 (eval = FALSE)
###################################################
##   panel.cumsum <- function(x,y,...){
##     x.ord <- order(x)
##     panel.xyplot(x[x.ord], cumsum(y[x.ord]), ...)}
## 
##    print(xyplot(Prob~NResp|factor(ClusterSize), groups=Trt, data=sh.mc, 
##        subset=ClusterSize>0&ClusterSize<13, type="s",
##        panel=panel.superpose, panel.groups=panel.cumsum,
##        as.table=T, auto.key=list(columns=4, lines=T, points=F),
## 			 xlab="Number of responses", ylab="Cumulative Probability R(R>=r|N=n)",
## 			 ylim=c(0,1.1)))


###################################################
### code chunk number 11: SOtest
###################################################
  set.seed(4461)
 (so.res <- trend.test(sh, test="SO", R=50, control=soControl(eps=0.1, max.directions=40)))


###################################################
### code chunk number 12: SOboot
###################################################
 hist(attr(so.res, "boot")$t[,1], freq=FALSE, xlab="Statistic", ylab="Density", main="")
 points(so.res$statistic, 0, pch="*",col="red", cex=3)


###################################################
### code chunk number 13: SOboot2 (eval = FALSE)
###################################################
##  hist(attr(so.res, "boot")$t[,1], freq=FALSE, xlab="Statistic", ylab="Density", main="")
##  points(so.res$statistic, 0, pch="*",col="red", cex=3)


###################################################
### code chunk number 14: RSGEEtest
###################################################
  trend.test(sh, test="RS")
  trend.test(sh, test="GEE")


###################################################
### code chunk number 15: SOmle
###################################################
   sh.SO.est <- SO.mc.est(sh, control=soControl(eps=0.1, max.directions=40))
   str(sh.SO.est)


###################################################
### code chunk number 16: SOmleplot
###################################################
   print(xyplot(Prob~NResp|factor(ClusterSize), groups=Trt, data=sh.SO.est, 
       subset=ClusterSize<13, type="s",
       panel=panel.superpose, panel.groups=panel.cumsum,
       as.table=T, auto.key=list(columns=4, lines=T, points=F),
			 xlab="Number of responses", ylab="Cumulative Probability R(R>=r|N=n)",
			 ylim=c(0,1.1), main=""))


###################################################
### code chunk number 17: SOmleplot2 (eval = FALSE)
###################################################
##    print(xyplot(Prob~NResp|factor(ClusterSize), groups=Trt, data=sh.SO.est, 
##        subset=ClusterSize<13, type="s",
##        panel=panel.superpose, panel.groups=panel.cumsum,
##        as.table=T, auto.key=list(columns=4, lines=T, points=F),
## 			 xlab="Number of responses", ylab="Cumulative Probability R(R>=r|N=n)",
## 			 ylim=c(0,1.1), main=""))


###################################################
### code chunk number 18: nostasot
###################################################
  NOSTASOT(sh, test="RS")


