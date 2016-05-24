## need a writable figs/ directory in getwd()
## -> switch to a temporary directory to save figures to
TMPDIR <- tempdir()
OWD <- setwd(TMPDIR)
dir.create("figs")


###################################################
### chunk number 1: 
###################################################
library("surveillance")
options(width=70)
options("prompt"="R> ")
set.seed(1234)

opendevice <- function(horizontal=TRUE,width=7,height=4,...) {
  #Do it for postscript instead  -- who uses postscript these days??
  args <- list(...)
  args$file <- sub(".pdf",".eps",args$file)
  args$width <- width
  args$height <- height
  args$horizontal <- FALSE
  do.call("postscript",args)

  par(mar=c(4,4,2,2))
}



###################################################
### chunk number 2: K1
###################################################
data("ha")
plot(aggregate(ha),main="Hepatitis A in Berlin 2001-2006")


###################################################
### chunk number 3: 
###################################################
opendevice(file="figs/002.pdf")
data("ha")
plot(aggregate(ha),main="Hepatitis A in Berlin 2001-2006")
dev.off()


###################################################
### chunk number 4: 
###################################################
sps <- sim.pointSource(p = 0.99, r = 0.5, length = 400, 
     A = 1, alpha = 1, beta = 0, phi = 0, frequency = 1, 
     state = NULL, K = 1.7)
plot(sps,xaxis.years=FALSE)


###################################################
### chunk number 5: 
###################################################
opendevice(file="figs/003.pdf")
plot(sps,xaxis.years=FALSE,legend.opts=list(x="topleft"))
dev.off()


###################################################
### chunk number 6: HAB662 eval=FALSE
###################################################
## ha.b662 <- algo.bayes(aggregate(ha), control = list(range = 209:290, b = 2, w = 6, alpha = 0.01))
## plot(ha.b662, firstweek=1, startyear = 2005)


###################################################
### chunk number 7: 
###################################################
ha.b662 <- algo.bayes(aggregate(ha), control = list(range = 209:290, b = 2, w = 6, alpha = 0.01))
plot(ha.b662, firstweek=1, startyear = 2005)
opendevice(file="figs/hab662.pdf")
plot(ha.b662, firstweek=1, startyear = 2005,legend.opts=list(x="topleft",horiz=TRUE))
dev.off()


###################################################
### chunk number 8: FACDC eval=FALSE
###################################################
## cntrl <- list(range = 300:400, m = 1, w = 3, b = 5, alpha = 0.01)
## sps.cdc <- algo.cdc(sps, control = cntrl)
## sps.farrington <- algo.farrington(sps, control = cntrl)


###################################################
### chunk number 9: 
###################################################
cntrl <- list(range = 300:400, m = 1, w = 3, b = 5, alpha = 0.01)
sps.cdc <- algo.cdc(sps, control = cntrl)
sps.farrington <- algo.farrington(sps, control = cntrl)


###################################################
### chunk number 10: 
###################################################
opendevice(file="figs/farringtoncdc.pdf")
par(mfcol = c(1, 2),cex=0.8)
plot(sps.cdc, legend = NULL, xaxis.years=FALSE)
plot(sps.farrington, legend = NULL, xaxis.years=FALSE)
dev.off()


###################################################
### chunk number 11: CUSUM eval=FALSE
###################################################
## kh <- find.kh(ARLa=500,ARLr=7)
## ha.cusum <- algo.cusum(aggregate(ha),control=list(k=kh$k,h=kh$h,m="glm",trans="rossi",range=209:290))


###################################################
### chunk number 12: 
###################################################
opendevice(file="figs/hacusum.pdf")
kh <- find.kh(ARLa=500,ARLr=7)
ha.cusum <- algo.cusum(aggregate(ha),control=list(k=kh$k,h=kh$h,m="glm",trans="rossi",range=209:290))
plot(ha.cusum,startyear=2005,legend.opts=list(x=30,y=5.5))
dev.off()

#Extract coefficients
beta <- coef(ha.cusum$control$m.glm)


###################################################
### chunk number 13: 
###################################################
print(algo.quality(ha.b662))


###################################################
### chunk number 14: 
###################################################
#This chunk contains stuff the reader should not see, but which is necessary 
#for the visual block to work.
control = list(
  list(funcName = "rki1"),
  list(funcName = "rki2"),
  list(funcName = "rki3"),
  list(funcName = "bayes1"),
  list(funcName = "bayes2"),
  list(funcName = "bayes3"),
#  list(funcName = "cdc",alpha=0.05,b=2,m=1),
#  list(funcName = "farrington",alpha=0.05,b=0,w=6),
  list(funcName = "farrington",alpha=0.05,b=1,w=6),
  list(funcName = "farrington",alpha=0.05,b=2,w=4))
control <- lapply(control,function(ctrl) {ctrl$range <- 300:400;return(ctrl)})

#Update range in each - cyclic continuation
data("k1")
range = (2*4*52) +  1:length(k1$observed)
aparv.control <- lapply(control,function(cntrl) { cntrl$range=range;return(cntrl)})

#Outbreaks
outbrks <- c("m1", "m2", "m3", "m4", "m5", "q1_nrwh", "q2", 
              "s1", "s2", "s3", "k1", "n1", "n2", "h1_nrwrp")

#Load and enlarge data.
outbrks <- lapply(outbrks,function(name) {
  #Load with data
  eval(substitute(data(name),list(name=name)))
  enlargeData(get(name),range=1:(4*52),times=2)
})

#Apply function to one
surv.one <- function(outbrk) {
  algo.compare(algo.call(outbrk,control=aparv.control))
}



###################################################
### chunk number 15:  eval=FALSE
###################################################
## #Apply function to one
## surv.one <- function(outbrk) {
##   algo.compare(algo.call(outbrk,control=aparv.control))
## }
## 
## algo.summary(lapply(outbrks, surv.one))
## 


###################################################
### chunk number 16: ALGOSUMMARY
###################################################
res <- algo.summary(lapply(outbrks,surv.one))


###################################################
### chunk number 17: 
###################################################
print(res,digits=3)


###################################################
### chunk number 18:  eval=FALSE
###################################################
## setClass( "sts", representation(week = "numeric",
##                                 freq = "numeric",
##                                 start = "numeric",
##                                 observed = "matrix",
##                                 state = "matrix",
##                                 alarm = "matrix",
##                                 upperbound  = "matrix",
##                                 neighbourhood= "matrix",
##                                 populationFrac= "matrix",
##                                 map = "SpatialPolygonsDataFrame",
##                                 control = "list"))
## 


###################################################
### chunk number 19: HA eval=FALSE
###################################################
## shp <- system.file("shapes/berlin.shp",package="surveillance")
## ha <- disProg2sts(ha, map=maptools::readShapePoly(shp,IDvar="SNAME"))
## plot(ha,type=observed ~ 1 | unit)
## 


###################################################
### chunk number 20: 
###################################################
opendevice(file="figs/ha-1unit.pdf",width=7,height=7)
par(mar=c(0,0,0,0))
shp <- system.file("shapes/berlin.shp",package="surveillance")
ha <- disProg2sts(ha, map=maptools::readShapePoly(shp,IDvar="SNAME"))
plot(ha,type=observed ~ 1 | unit)

dev.off()


###################################################
### chunk number 21: HA:MAP eval=FALSE
###################################################
## ha4 <- aggregate(ha[,c("pank","mitt","frkr","scho","chwi","neuk")],nfreq=13)
## ha4.cusum <- cusum(ha4,control=list(k=1.5,h=1.75,m="glm",trans="rossi",range=52:73))
## #ha4.b332 <- bayes(ha4,control=list(range=52:73,b=2,w=3,alpha=0.01/6))
## plot(ha4.cusum,type=observed ~ time | unit)


###################################################
### chunk number 22: 
###################################################
opendevice(file="figs/ha-timeunit.pdf",width=7,height=5)
ha4 <- aggregate(ha[,c("pank","mitt","frkr","scho","chwi","neuk")],nfreq=13)
ha4.cusum <- cusum(ha4,control=list(k=1.5,h=1.75,m="glm",trans="rossi",range=52:73))
#ha4.b332 <- bayes(ha4,control=list(range=52:73,b=2,w=3,alpha=0.01/6))
plot(ha4.cusum,type=observed ~ time | unit)
dev.off()


## finally switch back to original working directory
message("Note: selected figures have been saved in ", getwd(), "/figs")
setwd(OWD)
