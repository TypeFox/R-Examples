### R code from vignette source 'ThCointOverview.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: ThCointOverview.Rnw:37-39
###################################################
require(tsDyn)
options(prompt=" ", encoding="LATIN-9")


###################################################
### code chunk number 2: lib
###################################################
library(tsDyn)


###################################################
### code chunk number 3: grid1
###################################################
library(tsDyn)
data(lynx)
grid<-selectSETAR(lynx, m=1, thDelay=0, trim=0.15, criterion="SSR") 
print(grid)


###################################################
### code chunk number 4: plotgrid1
###################################################
plot(grid)


###################################################
### code chunk number 5: ThCointOverview.Rnw:549-551
###################################################
set<-setar(lynx, m=1, thDelay=0, th=grid$th)
summary(set)


###################################################
### code chunk number 6: grid2
###################################################
selectSETAR(lynx, m=1, thDelay=0, trim=0.15, criterion="SSR", nthresh=2) 


###################################################
### code chunk number 7: grid1
###################################################
selectSETAR(lynx, m=6, thDelay=0, trim=0.15, criterion="AIC", same.lags=TRUE)


###################################################
### code chunk number 8: ThCointOverview.Rnw:809-811 (eval = FALSE)
###################################################
## data(zeroyld)
## tvecm<-TVECM(zeroyld, nthresh=2,lag=1, ngridBeta=60, ngridTh=30, plot=TRUE,trim=0.05, beta=list(int=c(0.7, 1.1))) 


###################################################
### code chunk number 9: ThCointOverview.Rnw:969-971
###################################################
data(IIPUs)
set<-setar(IIPUs, m=16, thDelay=5, th=0.23)


###################################################
### code chunk number 10: ThCointOverview.Rnw:996-997 (eval = FALSE)
###################################################
## Hansen.Test<-setarTest(lynx, m=1, nboot=1000)


###################################################
### code chunk number 11: ThCointOverview.Rnw:1049-1054
###################################################
sun<-(sqrt(sunspot.year+1)-1)*2 
lin<-linear(sun, m=11)
set1<-setar(sun, m=11, th=7.4, thDelay=1, nested=TRUE)
set2<-setar(sun, m=11, th=c(5.3,8),nthresh=2, thDelay=1, nested=TRUE)
matrix(c(AIC(lin),AIC(set1),AIC(set2),BIC(lin),BIC(set1),BIC(set2)),ncol=2,dimnames=list(c("lin","set1", "set2"),c("AIC", "BIC")))


###################################################
### code chunk number 12: ThCointOverview.Rnw:1310-1314 (eval = FALSE)
###################################################
## data(zeroyld)
## dat<-zeroyld
## testSeo<-TVECM.SeoTest(dat, lag=1, beta=1, nboot=1000)
## summary(testSeo)


###################################################
### code chunk number 13: ThCointOverview.Rnw:1449-1454 (eval = FALSE)
###################################################
## system.time(test1<-TVECM.HStest(dat, lag=1, nboot=200))
## 
## library(doMC)
## registerDoMC(2) #Number of cores
## system.time(test1<-TVECM.HStest(dat, lag=1, nboot=200, hpc="foreach"))


