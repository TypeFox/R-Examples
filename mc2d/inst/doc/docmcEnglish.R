### R code from vignette source 'docmcEnglish.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: sh0
###################################################
set.seed(666)
options("width"=90,"digits"=3)


###################################################
### code chunk number 2: sh1
###################################################
library(mc2d)
ndvar(1001)
conc <- 10
cook <- mcstoc(rempiricalD, values=c(1,1/5,1/50), prob=c(0.027,0.373,0.600))
serving <- mcstoc(rgamma,shape=3.93,rate=0.0806)
expo <- conc * cook * serving
dose <- mcstoc(rpois,lambda=expo)
r <- 0.001
risk <- 1-(1-r)^dose
EC1 <- mc(cook,serving,expo,dose,risk)
print(EC1)
summary(EC1)


###################################################
### code chunk number 3: sh2
###################################################
ndunc(101)
conc <- mcstoc(rnorm,type="U",mean=10,sd=2)
cook <- mcstoc(rempiricalD, type="V",values=c(1,1/5,1/50), prob=c(0.027,0.373,0.600))
serving <- mcstoc(rgamma,type="V",shape=3.93,rate=0.0806)
expo <- conc * cook * serving
dose <- mcstoc(rpois,type="VU",lambda=expo)
r <- mcstoc(runif,type="U",min=0.0005,max=0.0015)
risk <- 1-(1-r)^dose
EC2 <- mc(conc,cook,serving,expo,dose,r,risk)
print(EC2,digits=2)
summary(EC2)


###################################################
### code chunk number 4: sh3 (eval = FALSE)
###################################################
## conc <- mcstoc(rnorm,type="U",mean=10,sd=2)
## cook <- mcstoc(rempiricalD, type="V",values=c(1,1/5,1/50), prob=c(0.027,0.373,0.600))
## serving <- mcstoc(rgamma,type="V",shape=3.93,rate=0.0806)
## ...
## dose <- mcstoc(rpois,type="VU",lambda=expo)
## r <- mcstoc(runif,type="U",min=0.0005,max=0.0015)
## ...


###################################################
### code chunk number 5: sh3bis
###################################################
x <- mcstoc(rnorm, mean=2, sd=3, rtrunc=TRUE, linf=1.5, lsup=2, lhs=TRUE)
summary(x)


###################################################
### code chunk number 6: sh3
###################################################
nu <- ndunc()
tmp <- (1:nu) > (nu/2)
mcdata(tmp,type="U")


###################################################
### code chunk number 7: sh4 (eval = FALSE)
###################################################
## ...
## expo <- conc * cook * serving
## ...
## risk <- 1-(1-r)^dose


###################################################
### code chunk number 8: sh4
###################################################
conc1 <- mcstoc(rnorm,type="U",mean=10,sd=2)
conc2 <- mcstoc(runif,type="U",min=8,max=12)
whichdist <- c(0.75,0.25)
concbis <- mcprobtree(whichdist,list("0"=conc1,"1"=conc2),type="U")


###################################################
### code chunk number 9: sh5
###################################################
cook < 1
suppressWarnings(tmp <- log(mcstoc(runif,min=-1,max=1)))
tmp
is.na(tmp)


###################################################
### code chunk number 10: sh5bis
###################################################
cornode(cook,serving,target=0.5,result=TRUE) 


###################################################
### code chunk number 11: sh6 (eval = FALSE)
###################################################
## ...
## EC2 <- mc(conc,cook,serving,expo,dose,r,risk)
## print(EC2)
## summary(EC2)


###################################################
### code chunk number 12: sh10
###################################################
modelEC3 <- mcmodel({
  conc <- mcstoc(rnorm,type="U",mean=10,sd=2)
  cook <- mcstoc(rempiricalD, type="V",values=c(1,1/5,1/50),
    prob=c(0.027,0.373,0.600))
  serving <- mcstoc(rgamma,type="V",shape=3.93,rate=0.0806)
  r <- mcstoc(runif,type="U",min=0.0005,max=0.0015)
  expo <- conc * cook * serving
  dose <- mcstoc(rpois,type="VU",lambda=expo)
  risk <- 1-(1-r)^dose
  mc(conc,cook,serving,expo,dose,r,risk)
})
modelEC3


###################################################
### code chunk number 13: sh11 (eval = FALSE)
###################################################
## EC3 <- evalmcmod(modelEC3,nsv=100,nsu=10,seed=666)
## EC4 <- evalmcmod(modelEC3,nsv=100,nsu=1000,seed=666)


###################################################
### code chunk number 14: sh12 (eval = FALSE)
###################################################
## modEC4 <- mcmodelcut({
## ## First block: unidimensional nodes
## {cook <- mcstoc(rempiricalD, type = "V", values = c(0, 1/5, 1/50), 
##                 prob = c(0.027, 0.373, 0.6))
##  serving <- mcstoc(rgamma, type = "V", shape = 3.93, rate = 0.0806)       
##  conc <- mcstoc(rnorm, type = "U", mean = 10, sd = 2)       
##  r <- mcstoc(runif, type = "U", min = 5e-04, max = 0.0015)     
## }
## ## Second block: two dimensional nodes
## {expo <- conc * cook * serving       
##  dose <- mcstoc(rpois, type = "VU", lambda = expo)
##  risk <- 1 - (1 - r)^dose
##  res <- mc(conc, cook, serving, expo, dose, r, risk)     }
## ## Third block: Outputs 
## {list(
##   sum = summary(res),
##   plot = plot(res, draw=FALSE),
##   minmax = lapply(res,range),
##   tor=tornado(res),
##   et =  sapply(res,sd))
## }
## })
## res <- evalmccut(modEC4, nsv = 10001, nsu = 101, seed = 666)
## summary(res)


###################################################
### code chunk number 15: sh14
###################################################
tmp <- summary(EC2,probs=c(0.995,0.999),digits=12)
tmp$risk


###################################################
### code chunk number 16: sh15 (eval = FALSE)
###################################################
## hist(EC2)


###################################################
### code chunk number 17: sh16
###################################################
hist(EC2)


###################################################
### code chunk number 18: sh17 (eval = FALSE)
###################################################
## plot(EC2)


###################################################
### code chunk number 19: sh18
###################################################
plot(EC2)


###################################################
### code chunk number 20: sh19
###################################################
torEC2 <- tornado(EC2)
plot(torEC2)


###################################################
### code chunk number 21: sh19b
###################################################
plot(torEC2)


###################################################
### code chunk number 22: sh19c
###################################################
tornadounc(EC2, output="risk", quant=.99)


###################################################
### code chunk number 23: sh19d
###################################################
mcratio(risk)


###################################################
### code chunk number 24: sh19ter
###################################################
tmp <- unmc(EC2, drop=TRUE)
dimu <- ncol(tmp$risk)
coef <- sapply(1:dimu, function(x) lm(tmp$risk[,x] ~ tmp$dose[,x])$coef) 
apply(coef,1,summary)


###################################################
### code chunk number 25: sh19ter
###################################################
mcstoc(runif, nvariates=3, min=c(1,2,3),max=4)


###################################################
### code chunk number 26: sh19quatr
###################################################
lim <- mcdata(c(1,2,3), type="0", nvariates=3)
mcstoc(runif, nvariates=3, min=lim,max=4)


###################################################
### code chunk number 27: sh20
###################################################
(p <- mcstoc(rdirichlet, type="U", nvariates=3, alpha=c(2,3,5)))
s <- mcstoc(rmultinomial,type="VU", nvariates=3, size=500, prob=p)
summary(s)


###################################################
### code chunk number 28: sh21
###################################################
sigma <- matrix(c(10,2,-5,2,10,-5,-5,-5,10), ncol=3)
(x <- mcstoc(rmultinormal,type="V", nvariates=3, mean=c(100,150,250), 
                              sigma=as.vector(sigma))) 
cor(x[,1,])


###################################################
### code chunk number 29: sh21
###################################################
m <- mcdata(c(100,150,250), type="0", nvariates=3)
mun <- mcstoc(rnorm, type="U", nvariates=3, mean=m, sd=20)
x <- mcstoc(rmultinormal, type="VU", nvariates=3, mean=mun, sigma=as.vector(sigma))
cor(x[,1,])


###################################################
### code chunk number 30: sh22
###################################################
val <- c(100,150,170,200)
pr <-  c(6,12,6,6)
out <- c('min','mean','max')
(x <- mcstoc(rempiricalD, type="U", outm=out, nvariates=30, 
              values=val,prob=pr))
mcstoc(rempiricalD,type="VU", values=x)


###################################################
### code chunk number 31: sh23
###################################################
conc1 <- mcstoc(rnorm, type="U", mean=10, sd=2)
conc2 <- mcstoc(runif, type="U", min=8, max=12)
conc <- mcdata(c(conc1,conc2),type="U",nvariates=2)

cook <- mcstoc(rempiricalD, type="V", values=c(1,1/5,1/50), prob=c(0.027,0.373,0.600))
serving <- mcstoc(rgamma,type="V",shape=3.93,rate=0.0806)
expo <- conc * cook * serving
dose <- mcstoc(rpois,type="VU",nvariates=2,lambda=expo)
r <- mcstoc(runif,type="U",min=0.0005,max=0.0015)
risk <- 1-(1-r)^dose
EC5 <- mc(conc,risk)
summary(EC5)


###################################################
### code chunk number 32: sh24
###################################################
mconc <- mcdata(c(10,14), type="0", nvariates=2)
conc <- mcstoc(rnorm, nvariates=2, type="U", mean=mconc, sd=2)

cook <- mcstoc(rempiricalD, type="V", values=c(1,1/5,1/50), prob=c(0.027,0.373,0.600))
serving <- mcstoc(rgamma,type="V",shape=3.93,rate=0.0806)
expo <- conc * cook * serving
dose <- mcstoc(rpois,type="VU",nvariates=2,lambda=expo)
dosetot <- mcapply(dose, margin="variates", fun=sum)
r <- mcstoc(runif,type="U",min=0.0005,max=0.0015)
risk <- 1-(1-r)^dosetot
EC6 <- mc(conc,risk)
summary(EC6)


###################################################
### code chunk number 33: Crypto1
###################################################
library(mc2d)
inca <- structure(c(0, 22.08, 60, 64.4, 72, 82.8, 90, 96, 100, 110, 120,  137.5, 144, 150, 160, 162.5, 165, 180, 182.5, 184, 192, 192.5,  200, 216, 220, 225, 230, 240, 250, 264, 270, 276, 288, 290, 300,  304, 312.8, 320, 322, 325, 330, 336, 340, 350, 360, 370, 375,  380, 384, 390, 400, 414, 420, 425, 430, 432, 432.5, 440, 442,  450, 456, 460, 460.8, 464, 470, 470.4, 480, 490, 500, 504, 510,  510.4, 516, 520, 525, 525.6, 528, 530, 540, 544, 550, 552, 560,  562, 565, 570, 576, 580, 582.5, 584, 585.6, 590, 596, 600, 606,  610, 614, 620, 625, 630, 635.4, 640, 648, 650, 656.2, 660, 664.4,  670, 670.4, 672, 675, 680, 682, 690, 696, 700, 710, 716, 720,  730, 730.4, 740, 744, 750, 756, 760, 774.8, 780, 784, 792, 796,  800, 810, 820, 828, 830, 840, 850, 850.4, 860, 864, 866.4, 870,  880, 890, 900, 908, 910, 915.2, 920, 930, 936, 950, 960, 970,  980, 984, 986.4, 990, 996, 1000, 1015.2, 1020, 1028, 1032, 1036,  1040, 1042, 1050, 1070, 1075, 1078.8, 1080, 1090, 1096, 1100,  1110, 1120, 1126.4, 1128, 1130, 1140, 1148, 1150, 1152, 1160,  1170, 1175, 1176.2, 1190, 1196, 1200, 1214, 1220, 1230, 1240,  1248, 1250, 1260, 1276, 1280, 1290, 1296, 1300, 1320, 1322, 1330,  1340, 1350, 1360, 1370, 1386.4, 1400, 1410, 1414, 1420, 1430,  1440, 1446, 1450, 1460, 1480, 1500, 1520, 1530, 1550, 1560, 1600,  1650, 1680, 1696, 1700, 1710, 1720, 1750, 1760, 1800, 1830, 1840,  1850, 1900, 1920, 1936, 1954, 1980, 1990, 2000, 2014, 2050, 2100,  2200, 2220, 2248, 2250, 2276, 2300, 2310, 2340, 2400, 2550, 2568,  2700, 2720, 2784, 2820, 2876, 3000, 3100, 3108, 3200, 2578, 7,  1, 8, 14, 3, 1, 1, 10, 1, 250, 1, 2, 120, 8, 6, 1, 5, 3, 12,  5, 5, 375, 2, 8, 7, 41, 408, 53, 4, 24, 7, 3, 2, 217, 1, 1, 44,  9, 1, 31, 1, 1, 17, 294, 5, 3, 9, 3, 12, 525, 5, 23, 1, 3, 4,  1, 28, 3, 154, 2, 5, 1, 2, 6, 1, 299, 8, 148, 1, 2, 1, 1, 8,  3, 1, 2, 14, 20, 1, 18, 2, 20, 6, 1, 8, 2, 8, 1, 1, 1, 4, 1,  487, 3, 5, 1, 7, 1, 5, 1, 24, 3, 17, 1, 42, 1, 2, 1, 1, 1, 16,  1, 3, 1, 30, 4, 1, 183, 4, 1, 5, 1, 141, 1, 14, 1, 12, 1, 2,  1, 206, 6, 2, 1, 4, 92, 10, 1, 5, 1, 3, 5, 5, 2, 87, 1, 1, 1,  5, 5, 4, 4, 78, 1, 3, 2, 1, 16, 1, 133, 1, 5, 1, 1, 1, 4, 1,  43, 1, 1, 1, 30, 1, 1, 7, 2, 6, 1, 1, 3, 3, 1, 10, 1, 5, 1, 1,  1, 1, 1, 159, 2, 1, 1, 10, 1, 16, 4, 2, 5, 3, 1, 3, 11, 4, 1,  2, 12, 5, 1, 1, 44, 3, 2, 1, 2, 17, 1, 4, 1, 1, 17, 1, 1, 3,  4, 18, 14, 4, 1, 2, 1, 1, 1, 2, 12, 1, 2, 1, 1, 1, 1, 1, 3, 1,  20, 1, 1, 1, 7, 1, 1, 3, 1, 2, 2, 1, 6, 1, 1, 1, 1, 1, 1, 1,  2, 1, 1, 2), .Dim = c(270L, 2L))
inca <- rep(inca[,1],inca[,2])/1000 


###################################################
### code chunk number 34: Crypto2
###################################################
hist(inca, xlab="l/day",freq=FALSE,main="") 


###################################################
### code chunk number 35: Crypto3
###################################################
ndvar(1001)
ndunc(101)
mcstoc(rempiricalD,type="V",values=inca)


###################################################
### code chunk number 36: Crypto4
###################################################
library(fitdistrplus) 
pzero <- sum(inca==0)/length(inca)
inca_non_0 <- inca[inca!=0] 
descdist(inca_non_0)


###################################################
### code chunk number 37: Crypto5
###################################################
descdist(inca_non_0)


###################################################
### code chunk number 38: Crypto6
###################################################
Adj_water <- fitdist(inca_non_0,"lnorm",method="mle")
meanlog <- Adj_water$est[1]
sdlog   <- Adj_water$est[2]
summary(Adj_water)


###################################################
### code chunk number 39: Crypto7bis (eval = FALSE)
###################################################
## Boot <- bootdist(Adj_water, bootmethod="param", niter=ndunc())
## Mean_conso <- mcdata(Boot$estim$meanlog, type="U")
## Sd_conso <- mcdata(Boot$estim$sdlog, type="U")
## conso1 <- mcstoc(rlnorm, type="VU", meanlog= Mean_conso, sdlog= Sd_conso)


###################################################
### code chunk number 40: Crypto8
###################################################
conso0 <- mcdata(0,type="V")
conso1 <- mcstoc(rlnorm, type="V", meanlog=meanlog, sdlog=sdlog)
v <- mcprobtree(c(pzero,1-pzero), list("0"=conso0,"1"=conso1), type = "V")
summary(v)


###################################################
### code chunk number 41: Crypto9
###################################################
datDR <- list(   dose=c(30,100,300,500,1000,10000,100000,1000000),
                 pi=c(2,4,2,5,2,3,1,1),
                 ni=c(5,8,3,6,2,3,1,1))

estDR <- function(pos,ref){
   suppressWarnings(
	-glm(cbind(ref$ni-pos,pos) ~ ref$dose + 0, 
		binomial(link="log"))$coefficients)}

ml <- 1-exp(-estDR(datDR$pi, datDR) * datDR$dose)

DR <- function(n){
   boot <- matrix(rbinom(length(datDR$dose)*n,datDR$ni,ml),nrow=length(datDR$dose))
   apply(boot,2,estDR,ref=datDR)}
r <- mcstoc(DR, type="U")
summary(r)


###################################################
### code chunk number 42: Crypto10
###################################################
Rr <- mcstoc(rbeta, type="U", shape1=2.65, shape2=3.64)
w <- mcstoc(rbeta, type="V", shape1=2.6, shape2=3.4) 


###################################################
### code chunk number 43: Crypto11
###################################################
Oo <- 2
l <- (Oo + mcstoc(rnbinom, type="U", size=Oo+1, prob=Rr))/100


###################################################
### code chunk number 44: Crypto12
###################################################
Or <- l * v * w
P <- 10000 * (1-exp(-r*Or))
summary(P)


###################################################
### code chunk number 45: Crypto13 (eval = FALSE)
###################################################
## Oo <- mcdata(c(0,1,2,5,10,20,50,100,1000),type="0",nvariates=9) 


###################################################
### code chunk number 46: Crypto11
###################################################
Oo <- 2
l <- (Oo + mcstoc(rnbinom, type="U", size=Oo+1, prob=Rr))/100


###################################################
### code chunk number 47: Crypto12
###################################################
Or <- l * v * w
P <- 10000 * (1-exp(-r*Or))
summary(P)


###################################################
### code chunk number 48: Crypto13 (eval = FALSE)
###################################################
## Oo <- mcdata(c(0,1,2,5,10,20,50,100,1000),type="0",nvariates=9) 


