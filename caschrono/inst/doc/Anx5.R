### R code from vignette source 'Anx5.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Anx5.Rnw:136-140
###################################################
owidth <- getOption("width") # largeur des sorties
options(width=60, continue="+ ","warn"=-1 )
.PngNo <- 0
nom.fich = "anx5-bitmap-"


###################################################
### code chunk number 2: bfig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 7, height = 7, pointsize = 12, bg = "white")


###################################################
### code chunk number 3: bfigps (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""), width = 7, height = 7, pointsize = 12, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 4: bfig1 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 5, height = 2, pointsize = 10, bg = "white")


###################################################
### code chunk number 5: bfigps1 (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""),  width = 5, height =2, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 6: bfig2 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 3.9, height = 3.1, pointsize = 10, bg = "white")


###################################################
### code chunk number 7: bfigps2 (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""), width = 3.9, height = 3.1,   pointsize = 10, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 8: bfig3 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 5.92, height = 6.74, pointsize = 10, bg = "white")


###################################################
### code chunk number 9: bfigps3 (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""), width = 5.92, height = 6.74, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 10: bfig4 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 6, height = 6, pointsize = 10, bg = "white")


###################################################
### code chunk number 11: bfigps4 (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""), width = 6, height = 6, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 12: zfig2 (eval = FALSE)
###################################################
## dev.null <- dev.off()


###################################################
### code chunk number 13: zfiginclude (eval = FALSE)
###################################################
## cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 14: period
###################################################
# serie de 48 points
# x1 et x2 sont les fonctions periodiques
# et on verifie que leurs differences saisonnieres sont nulles.
temps = 1:48
x1 = cos(2*pi*temps/4)
diff(x1,4)
x2 = sin(2*pi*temps/4)
diff(x2,4)


###################################################
### code chunk number 15: Anx5.Rnw:267-277
###################################################
require(dse)
require(polynom)
require(forecast)
nobs=200
set.seed(234)
ar.1=c(1,.8)
masaiso = polynomial(c(1,rep(0,11),.5))
ar.13=  polynomial(c(1,rep(0,11),-1))*polynomial(ar.1)
ar.2 = polynomial(c(1,-1))*polynomial(ar.1)
moy=5


###################################################
### code chunk number 16: arimasim
###################################################
ya = arima.sim(n=nobs,list(order=c(1,0,12),ar=-.8,
ma=c(rep(0,11),.5))) + moy
mean(ya)
Arima(ya,order=c(1,0,0),seasonal=list(order=c(0,0,1),period=12),
include.mean=TRUE)


###################################################
### code chunk number 17: arimasim2
###################################################
# I(1)
yd0=diffinv(ya,1)[-1]
x1 = as.matrix(seq(1,length(yd0)))
Arima(yd0,order=c(1,1,0),
seasonal=list(order=c(0,0,1), period=12),xreg=x1)
Arima(yd0,order=c(1,1,0),
seasonal=list(order=c(0,0,1), period=12),
include.drift=TRUE)


###################################################
### code chunk number 18: arimasim3
###################################################
# I(12)
yd2=diffinv(yd0,12)[-(1:12)]
x2 = as.matrix(seq(1,length(yd0))^2)
(m1=Arima(yd2,order=c(1,1,0),seasonal=list(order=c(0,1,1),
         frequency=12),xreg=x2))
(m1b=Arima(yd2,order=c(1,1,0),seasonal=list(order=c(0,1,1),
       frequency=12),include.drift=TRUE))


###################################################
### code chunk number 19: simulate1
###################################################
(ar.1=polynomial(c(1,.8)))
(ar.2 = polynomial(c(1,-1))*ar.1)
# terme autoregressif complet jusqu'au retard 14
(ar.14 =  polynomial(c(1,rep(0,11),-1))*ar.2)
moy=5   
(masaiso = polynomial(c(1,rep(0,11),.5)))
(MA2= array(masaiso,c(13,1,1)))
# 
AR2 = array(ar.14,c(15,1,1))
MA2= array(masaiso,c(13,1,1))
modsarima=   ARMA(A=AR2, B=MA2,C=1)
derive= moy
cte =  derive*predict(polynomial(ar.1),1)
u.t = cte*matrix(rep(1,nobs))
y2 = simulate(modsarima,input=u.t,sampleT=nobs)$output


###################################################
### code chunk number 20: Anx5.Rnw:363-364
###################################################
mean(diff(diff(y2,1),12))


###################################################
### code chunk number 21: Anx5.Rnw:370-381
###################################################
x1 = as.matrix(seq(1,length(y2)))
x2 = as.matrix(seq(1,length(y2))^2)
(m1=Arima(y2,order=c(1,1,0),
seasonal=list(order=c(0,1,1), period=12), 
xreg=x1))
(m2=Arima(y2,order=c(1,1,0),
seasonal=list(order=c(0,1,1),period=12),
include.drift=TRUE))
(m3=Arima(y2,order=c(1,1,0),
seasonal=list(order=c(0,1,1),period=12),
xreg=x2))


###################################################
### code chunk number 22: Anx5.Rnw:385-386
###################################################
cf3 = m3$coef


###################################################
### code chunk number 23: derive (eval = FALSE)
###################################################
## set.seed(51)
## c=-.2
## y0=arima.sim(list(order=c(1,0,0),ar=.9),n=200)+c
## y1=diffinv(y0,4)[-(1:4)]
## lag.plot(rev(y1),4,layout=c(2,2),do.lines=FALSE,
##    main="Derive negative",diag.col="red")
## c=.2
## y0=arima.sim(list(order=c(1,0,0),ar=.9),n=200)+c
## y1=diffinv(y0,4)[-(1:4)]
## lag.plot(rev(y1),4,layout=c(2,2),ask=NA,do.lines=FALSE,
##                      main="Derive positive",diag.col="red")


###################################################
### code chunk number 24: spur1
###################################################
set.seed(514) ; nobs=300
y0 = 2+rnorm(nobs)
y = diffinv(y0)
x0 = 2-1.5*rnorm(nobs)
x = diffinv(x0)


###################################################
### code chunk number 25: spur20 (eval = FALSE)
###################################################
## plot(x,y,type=l)


###################################################
### code chunk number 26: spur2 (eval = FALSE)
###################################################
## plot.ts(cbind(x,y), plot.type = "single",lty=1:2)
## legend("topleft",c("x","y"),lty=1:2)


###################################################
### code chunk number 27: spur1c (eval = FALSE)
###################################################
## op <- par(mfrow=c(2,1), mar=c(4,4,1,1), oma=c(0,0,0,0))
## plot(x,y, xlab='x', ylab='y',pch='+',cex=.6)
## plot.ts(cbind(x,y), plot.type='single', ylab='x, y', xlab='temps')
## par(op)


###################################################
### code chunk number 28: spur2.a
###################################################
.PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
pdf(file=paste(file,".pdf",sep=""), width = 3.9, height = 3.1, pointsize = 10, bg = "white")
op <- par(mfrow=c(2,1), mar=c(4,4,1,1), oma=c(0,0,0,0))
plot(x,y, xlab='x', ylab='y',pch='+',cex=.6)
plot.ts(cbind(x,y), plot.type='single', ylab='x, y', xlab='temps')
par(op)
dev.null <- dev.off()
postscript(file=paste(file,".ps",sep=""), width = 3.9, height = 3.1,   pointsize = 10, bg = "white",horizontal= FALSE,paper="special")
op <- par(mfrow=c(2,1), mar=c(4,4,1,1), oma=c(0,0,0,0))
plot(x,y, xlab='x', ylab='y',pch='+',cex=.6)
plot.ts(cbind(x,y), plot.type='single', ylab='x, y', xlab='temps')
par(op)
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 29: spur1
###################################################
mod4 = lm(y~x)
(aa=summary(mod4))


###################################################
### code chunk number 30: spur4 (eval = FALSE)
###################################################
## acf(resid)


###################################################
### code chunk number 31: spur2d (eval = FALSE)
###################################################
## op <- par(mfrow=c(2,1), mar=c(4,3,1,0), oma=c(0,0,0,0))
## plot.ts(mod4$residuals, xlab='temps',ylab='residu MCO')
## abline(h=0)
## acf(mod4$residuals,main="")
## par(op)


###################################################
### code chunk number 32: spur2.b
###################################################
.PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
pdf(file=paste(file,".pdf",sep=""), width = 7, height = 7, pointsize = 12, bg = "white")
op <- par(mfrow=c(2,1), mar=c(4,3,1,0), oma=c(0,0,0,0))
plot.ts(mod4$residuals, xlab='temps',ylab='residu MCO')
abline(h=0)
acf(mod4$residuals,main="")
par(op)
dev.null <- dev.off()
postscript(file=paste(file,".ps",sep=""), width = 7, height = 7, pointsize = 12, bg = "white",horizontal= FALSE,paper="special")
op <- par(mfrow=c(2,1), mar=c(4,3,1,0), oma=c(0,0,0,0))
plot.ts(mod4$residuals, xlab='temps',ylab='residu MCO')
abline(h=0)
acf(mod4$residuals,main="")
par(op)
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


