### R code from vignette source 'kedd.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: kedd.Rnw:132-133
###################################################
library(kedd)


###################################################
### code chunk number 2: kedd.Rnw:135-137
###################################################
kernel.fun(x = seq(-0.02,0.02,by=0.01), deriv.order = 1, kernel = "gaussian")$kx
kernel.conv(x = seq(-0.02,0.02,by=0.01), deriv.order = 1, kernel = "gaussian")$kx


###################################################
### code chunk number 3: kedd.Rnw:139-141
###################################################
plot(kernel.fun(deriv.order = 1, kernel = "gaussian"))
plot(kernel.conv(deriv.order = 1, kernel = "gaussian"))


###################################################
### code chunk number 4: kedd.Rnw:146-147
###################################################
plot(kernel.fun(deriv.order = 1, kernel = "gaussian"),main = "", sub = "")


###################################################
### code chunk number 5: kedd.Rnw:149-150
###################################################
plot(kernel.conv(deriv.order = 1, kernel = "gaussian"),main = "", sub = "")


###################################################
### code chunk number 6: kedd.Rnw:176-182
###################################################
fx <- function(x) 0.5 *(-4*x-6)* dnorm(x,-1.5,0.5) + 0.5 *(-4*x+6) * dnorm(x,1.5,0.5)
kernels <- c("gaussian","biweight","triweight","tricube")
plot(dkde(x = bimodal, deriv.order = 1, h = 0.6, kernel = kernels[1]),col = 1,ylim=c(-0.6,0.6) ,sub="", main="")
for (i in 2:length(kernels))lines(dkde(x = bimodal, deriv.order = 1, h = 0.6, kernel = kernels[i] ), col = i)
curve(fx,add=TRUE,lty=8)
legend("topright", legend = c(TRUE,kernels), col = c("black",seq(kernels)),lty = c(8,rep(1,length(kernels))), inset = .015)


###################################################
### code chunk number 7: kedd.Rnw:184-189
###################################################
h <- c(0.14,0.3,0.6,1.2)                                                                                                               
plot(dkde(x = bimodal, deriv.order = 1, h = h[1], kernel = kernels[1]),col = 1,ylim=c(-0.6,1) ,sub="", main="")                   
for (i in 2:length(h))lines(dkde(x = bimodal, deriv.order = 1, h = h[i], kernel = kernels[1] ), col = i)                            
curve(fx,add=TRUE,lty=8)                                                                                                            
legend("topright", legend = c("TRUE",paste("h =",bquote(.(h)))), col = c("black",seq(h)),lty = c(8,rep(1,length(h))), inset = .015) 


###################################################
### code chunk number 8: kedd.Rnw:229-233
###################################################
hatf  <- dkde(bimodal, deriv.order = 0)
hatf1 <- dkde(bimodal, deriv.order = 1)
hatf2 <- dkde(bimodal, deriv.order = 2)
hatf3 <- dkde(bimodal, deriv.order = 3)


###################################################
### code chunk number 9: kedd.Rnw:237-248
###################################################
fx  <- function(x) 0.5 * dnorm(x,-1.5,0.5) + 0.5 * dnorm(x,1.5,0.5)
fx1 <- function(x) 0.5 *(-4*x-6)* dnorm(x,-1.5,0.5) + 0.5 *(-4*x+6) * 
                   dnorm(x,1.5,0.5)
fx2 <- function(x) 0.5 * ((-4*x-6)^2 - 4) * dnorm(x,-1.5,0.5) + 0.5 *
                   ((-4*x+6)^2 - 4) * dnorm(x,1.5,0.5)
fx3 <- function(x) 0.5 * (-4*x-6) * ((-4*x-6)^2 - 12) * dnorm(x,-1.5,0.5) +
                     0.5 * (-4*x+6) * ((-4*x+6)^2 - 12) * dnorm(x,1.5,0.5)
plot(hatf ,fx = fx)
plot(hatf1,fx = fx1)
plot(hatf2,fx = fx2)
plot(hatf3,fx = fx3)


###################################################
### code chunk number 10: kedd.Rnw:253-254
###################################################
plot(hatf,fx = fx,lwd=2,sub="",main="")


###################################################
### code chunk number 11: kedd.Rnw:256-257
###################################################
plot(hatf1,fx = fx1,lwd=2,sub="",main="")


###################################################
### code chunk number 12: kedd.Rnw:259-260
###################################################
plot(hatf2,fx = fx2,lwd=2,sub="",main="")


###################################################
### code chunk number 13: kedd.Rnw:262-263
###################################################
plot(hatf3,fx = fx3,lwd=2,sub="",main="")


###################################################
### code chunk number 14: kedd.Rnw:312-316
###################################################
h.amise(bimodal, deriv.order = 0)
h.amise(bimodal, deriv.order = 1)
h.amise(bimodal, deriv.order = 2)
h.amise(bimodal, deriv.order = 3)


###################################################
### code chunk number 15: kedd.Rnw:352-356
###################################################
kernels <- eval(formals(h.mlcv.default)$kernel)
hmlcv <- numeric()
for(i in 1:length(kernels))  
             hmlcv[i] <- h.mlcv(bimodal, kernel =  kernels[i])$h


###################################################
### code chunk number 16: kedd.Rnw:358-359
###################################################
data.frame(kernels,hmlcv)


###################################################
### code chunk number 17: kedd.Rnw:362-364
###################################################
plot(h.mlcv(bimodal, kernel =  kernels[1]), seq.bws = seq(0.1,1,length=50))
plot(h.mlcv(bimodal, kernel =  kernels[2]), seq.bws = seq(0.1,1,length=50))


###################################################
### code chunk number 18: kedd.Rnw:369-370
###################################################
plot(h.mlcv(bimodal, kernel =  kernels[1]), seq.bws = seq(0.1,1,length=50),sub="",main="")


###################################################
### code chunk number 19: kedd.Rnw:372-373
###################################################
plot(h.mlcv(bimodal, kernel =  kernels[2]), seq.bws = seq(0.1,1,length=50),sub="",main="")


###################################################
### code chunk number 20: kedd.Rnw:412-416
###################################################
h.ucv(bimodal, deriv.order = 0)
h.ucv(bimodal, deriv.order = 1)
h.ucv(bimodal, deriv.order = 2)
h.ucv(bimodal, deriv.order = 3)


###################################################
### code chunk number 21: kedd.Rnw:419-420
###################################################
for (i in 0:3) plot(h.ucv(bimodal, deriv.order = i))


###################################################
### code chunk number 22: kedd.Rnw:425-426
###################################################
plot(h.ucv(bimodal, deriv.order = 0),sub="",main="")


###################################################
### code chunk number 23: kedd.Rnw:428-429
###################################################
plot(h.ucv(bimodal, deriv.order = 1),sub="",main="")


###################################################
### code chunk number 24: kedd.Rnw:431-432
###################################################
plot(h.ucv(bimodal, deriv.order = 2),sub="",main="")


###################################################
### code chunk number 25: kedd.Rnw:434-435
###################################################
plot(h.ucv(bimodal, deriv.order = 3),sub="",main="")


###################################################
### code chunk number 26: kedd.Rnw:481-485
###################################################
h.bcv(bimodal, whichbcv = 1, deriv.order = 0)
h.bcv(bimodal, whichbcv = 2, deriv.order = 0)
h.bcv(bimodal, whichbcv = 1, deriv.order = 1, lower=0.1, upper=0.8)
h.bcv(bimodal, whichbcv = 2, deriv.order = 1, lower=0.1, upper=0.8)


###################################################
### code chunk number 27: kedd.Rnw:488-499
###################################################
## deriv.order = 0
plot(h.bcv(bimodal, whichbcv = 2, deriv.order = 0))
lines(h.bcv(bimodal, whichbcv = 1, deriv.order = 0),col="red")
legend("topright", c("BCV1","BCV2"),lty=1,col=c("red","black"),
       inset = .015)
## deriv.order = 1
plot(h.bcv(bimodal, whichbcv = 2, deriv.order = 1),seq.bws = 
     seq(0.1,0.8,length=50))
lines(h.bcv(bimodal, whichbcv = 1, deriv.order = 1),col="red")
legend("topright", c("BCV1","BCV2"),lty=1,col=c("red","black"),
      inset = .015)


###################################################
### code chunk number 28: kedd.Rnw:504-507
###################################################
plot(h.bcv(bimodal, whichbcv = 2, deriv.order = 0),sub="",main="")
lines(h.bcv(bimodal, whichbcv = 1, deriv.order = 0),col="red")
legend("topright", c("BCV1","BCV2"),lty=1,col=c("red","black"),inset = .015)


###################################################
### code chunk number 29: kedd.Rnw:509-512
###################################################
plot(h.bcv(bimodal, whichbcv = 2, deriv.order = 1),seq.bws=seq(0.1,0.8,length=50),sub="",main="")
lines(h.bcv(bimodal, whichbcv = 1, deriv.order = 1),col="red")
legend("topright", c("BCV1","BCV2"),lty=1,col=c("red","black"),inset = .015)


###################################################
### code chunk number 30: kedd.Rnw:553-557
###################################################
h.ccv(bimodal, deriv.order = 0, upper = 0.5)
h.ccv(bimodal, deriv.order = 1, upper = 0.5)
h.ccv(bimodal, deriv.order = 2, upper = 0.5)
h.ccv(bimodal, deriv.order = 3, upper = 0.5)


###################################################
### code chunk number 31: kedd.Rnw:560-562
###################################################
for (i in 0:3) 
 plot(h.ccv(bimodal, deriv.order = i), seq.bws=seq(0.1,0.5,length=50))


###################################################
### code chunk number 32: kedd.Rnw:567-568
###################################################
plot(h.ccv(bimodal, deriv.order = 0),seq.bws=seq(0.1,0.5,length=50),sub="",main="")


###################################################
### code chunk number 33: kedd.Rnw:570-571
###################################################
plot(h.ccv(bimodal, deriv.order = 1),seq.bws=seq(0.1,0.5,length=50),sub="",main="")


###################################################
### code chunk number 34: kedd.Rnw:573-574
###################################################
plot(h.ccv(bimodal, deriv.order = 2),seq.bws=seq(0.1,0.5,length=50),sub="",main="")


###################################################
### code chunk number 35: kedd.Rnw:576-577
###################################################
plot(h.ccv(bimodal, deriv.order = 3),seq.bws=seq(0.1,0.5,length=50),sub="",main="")


###################################################
### code chunk number 36: kedd.Rnw:614-618
###################################################
h.mcv(bimodal, deriv.order = 0, upper = 0.5)
h.mcv(bimodal, deriv.order = 1, upper = 0.5)
h.mcv(bimodal, deriv.order = 2, upper = 0.5)
h.mcv(bimodal, deriv.order = 3, upper = 0.5)


###################################################
### code chunk number 37: kedd.Rnw:621-623
###################################################
for (i in 0:3)
 plot(h.mcv(bimodal, deriv.order = i), seq.bws=seq(0.1,0.5,length=50))


###################################################
### code chunk number 38: kedd.Rnw:628-629
###################################################
plot(h.mcv(bimodal, deriv.order = 0),seq.bws=seq(0.1,0.5,length=50),sub="",main="")


###################################################
### code chunk number 39: kedd.Rnw:631-632
###################################################
plot(h.mcv(bimodal, deriv.order = 1),seq.bws=seq(0.1,0.5,length=50),sub="",main="")


###################################################
### code chunk number 40: kedd.Rnw:634-635
###################################################
plot(h.mcv(bimodal, deriv.order = 2),seq.bws=seq(0.1,0.5,length=50),sub="",main="")


###################################################
### code chunk number 41: kedd.Rnw:637-638
###################################################
plot(h.mcv(bimodal, deriv.order = 3),seq.bws=seq(0.1,0.5,length=50),sub="",main="")


###################################################
### code chunk number 42: kedd.Rnw:676-680
###################################################
h.tcv(bimodal, deriv.order = 0)
h.tcv(bimodal, deriv.order = 1)
h.tcv(bimodal, deriv.order = 2)
h.tcv(bimodal, deriv.order = 3)


###################################################
### code chunk number 43: kedd.Rnw:683-684
###################################################
for (i in 0:3) plot(h.tcv(bimodal, deriv.order = i))


###################################################
### code chunk number 44: kedd.Rnw:689-690
###################################################
plot(h.tcv(bimodal, deriv.order = 0),seq.bws=seq(0.1,1.5,length=50),sub="",main="")


###################################################
### code chunk number 45: kedd.Rnw:692-693
###################################################
plot(h.tcv(bimodal, deriv.order = 1),seq.bws=seq(0.3,1.5,length=50),sub="",main="")


###################################################
### code chunk number 46: kedd.Rnw:695-696
###################################################
plot(h.tcv(bimodal, deriv.order = 2),seq.bws=seq(0.5,1.5,length=50),sub="",main="")


###################################################
### code chunk number 47: kedd.Rnw:698-699
###################################################
plot(h.tcv(bimodal, deriv.order = 3),seq.bws=seq(0.5,1.5,length=50),sub="",main="")


