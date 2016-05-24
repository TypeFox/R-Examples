### R code from vignette source 'Anx4.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Anx4.Rnw:138-142
###################################################
owidth <- getOption("width") # largeur des sorties
options(width=60, continue="+ ","warn"=-1 )
.PngNo <- 0
nom.fich = "anx4-bitmap-"


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
### code chunk number 14: lb2
###################################################
require(fBasics)
require(caschrono)
data(csdl)
aa = returns(csdl, percentage = TRUE)
aab = aa[complete.cases(aa) == TRUE,]
r.csdl = its(aab, as.POSIXct(row.names(aab)))
r.danone= r.csdl[,3]
rdt2 = (r.danone-mean(r.danone))^2
a0=Box.test.2(r.danone, nlag=c(3,6,9,12),
 type="Ljung-Box",decim=4)
a1=Box.test.2(r.danone[1:600], nlag=c(3,6,9,12),
 type="Ljung-Box",decim=4)
a2=Box.test.2(rdt2[1:600], nlag=c(3,6,9,12), 
 type="Ljung-Box",decim=4)
a12 = cbind(a0,a1[,2],a2[,2])
colnames(a12)= c("Retard","p-val. serie compl.",
 "p-val. 600 obs.", "p-val. rdt carre")


###################################################
### code chunk number 15: lb2.tex
###################################################
require(xtable)
xtable(a12, caption="Table a12 : test de blancheur 
du rendement de Danone et de son carre.",label="danoblanc", 
digits=4)


###################################################
### code chunk number 16: Anx4.Rnw:284-286
###################################################
ARMAtoMA(-.7, 0, 10)
ARMAtoMA(c(rep(0,11),-.7), 0, 25)


###################################################
### code chunk number 17: Anx4.Rnw:293-296
###################################################
require(FitARMA)
ImpulseCoefficientsARMA(-.7,0,lag.max=10)
ImpulseCoefficientsARMA(c(rep(0,11),-.7),0,lag.max=25)


###################################################
### code chunk number 18: prep_ma2
###################################################
set.seed(951)
ya = arima.sim(n = 200, list( ma = c(-0.3, 0.6)), sd = sqrt(1.5))


###################################################
### code chunk number 19: plot.ma2 (eval = FALSE)
###################################################
## set.seed(951)
## ya = arima.sim(n = 200, list( ma = c(-0.3, 0.6)), sd = sqrt(1.5))
## titre= "MA(2)"
## plotacfthemp(ya, ma=c(-0.3, 0.6), lag.max=20)


###################################################
### code chunk number 20: Anx4.Rnw:352-359
###################################################
.PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
pdf(file=paste(file,".pdf",sep=""), width = 6, height = 6, pointsize = 10, bg = "white")
set.seed(951)
ya = arima.sim(n = 200, list( ma = c(-0.3, 0.6)), sd = sqrt(1.5))
titre= "MA(2)"
plotacfthemp(ya, ma=c(-0.3, 0.6), lag.max=20)
dev.null <- dev.off()
postscript(file=paste(file,".ps",sep=""), width = 6, height = 6, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")
set.seed(951)
ya = arima.sim(n = 200, list( ma = c(-0.3, 0.6)), sd = sqrt(1.5))
titre= "MA(2)"
plotacfthemp(ya, ma=c(-0.3, 0.6), lag.max=20)
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 21: ma2.estim
###################################################
(mod.ma2= Arima(ya,order=c(0,0,2), include.mean=FALSE))


###################################################
### code chunk number 22: prep_ar1
###################################################
set.seed(5419)
n2 = 210
yb = arima.sim(n=200, list(ar=-0.8), sd=sqrt(1.5))
yb = yb-10


###################################################
### code chunk number 23: plot.ar1 (eval = FALSE)
###################################################
## set.seed(5419)
## n2 = 210
## yb = arima.sim(n=200, list(ar=-0.8), sd=sqrt(1.5))
## yb = yb-10
## plotacfthemp(yb, ar=-0.8, lag.max=20)


###################################################
### code chunk number 24: Anx4.Rnw:396-403
###################################################
.PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
pdf(file=paste(file,".pdf",sep=""), width = 7, height = 7, pointsize = 12, bg = "white")
set.seed(5419)
n2 = 210
yb = arima.sim(n=200, list(ar=-0.8), sd=sqrt(1.5))
yb = yb-10
plotacfthemp(yb, ar=-0.8, lag.max=20)
dev.null <- dev.off()
postscript(file=paste(file,".ps",sep=""), width = 7, height = 7, pointsize = 12, bg = "white",horizontal= FALSE,paper="special")
set.seed(5419)
n2 = 210
yb = arima.sim(n=200, list(ar=-0.8), sd=sqrt(1.5))
yb = yb-10
plotacfthemp(yb, ar=-0.8, lag.max=20)
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 25: estim.12
###################################################
(mod12 = Arima(yb, order=c(1,0,0)))


###################################################
### code chunk number 26: simar2
###################################################
require(FitARMA)
set.seed(51)
nsim=100
nobs=200
nsim=50
nlag=20
y.mat = matrix(0,nrow=nobs,ncol=nsim)
facp.mat= matrix(0,nrow=nlag,ncol=nsim)
y.mat = matrix(0,nrow=200,ncol=nsim)
facp.mat= matrix(0,nrow=nlag,ncol=nsim)
for (isim in 1:nsim)
{
# isim=2
y.mat[,isim] = arima.sim(n=nobs,list(ar=c(-0.7, 0.2)),sd = sqrt(2))
facp.mat[,isim] = pacf(y.mat[,isim], 20, plot=FALSE)$acf
}
aa=t(apply(facp.mat,1,'quantile',probs = c(0.25,.75)))
#pacf theo
theo =  TacvfARMA(phi=c(-.7,.2),lag.max=20)
pacf.theo  = PacfDL(theo/theo[1], LinearPredictor=TRUE)$Pacf
# intervalle autour de 0 \`a 50%
binf= qnorm(.25)/nobs^.5
bsup=qnorm(.75)/nobs^.5
aaa = cbind(aa,pacf.theo,binf,bsup)


###################################################
### code chunk number 27: Anx4.Rnw:483-485 (eval = FALSE)
###################################################
## matplot(1:20,aaa,type='l',ylab="PACF",xlab="retard",col='black')
## legend("topright",paste("nombre de simulations : ",as.character(nsim)))


###################################################
### code chunk number 28: matplot0 (eval = FALSE)
###################################################
## matplot(1:20,aaa,type='l',ylab="PACF",xlab="retard",col='black')
## legend("topright",paste("nombre de simulations : ",as.character(nsim)))


###################################################
### code chunk number 29: Anx4.Rnw:500-507
###################################################
.PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
pdf(file=paste(file,".pdf",sep=""), width = 6, height = 6, pointsize = 10, bg = "white")
matplot(1:20,aaa,type='l',ylab="PACF",xlab="retard",col='black')
legend("topright",paste("nombre de simulations : ",as.character(nsim)))
dev.null <- dev.off()
postscript(file=paste(file,".ps",sep=""), width = 6, height = 6, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")
matplot(1:20,aaa,type='l',ylab="PACF",xlab="retard",col='black')
legend("topright",paste("nombre de simulations : ",as.character(nsim)))
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 30: prep.y1
###################################################
set.seed(123)
y1 = arima.sim(n=100,list(ar=-.7), sd=sqrt(4))
y2 = arima.sim(n=100,list(ar=c(rep(0,11),-.7)), sd=sqrt(4))


###################################################
### code chunk number 31: estim.ma1
###################################################
(mod2=Arima(y1,order=c(0,0,1), include.mean=FALSE))


###################################################
### code chunk number 32: test.ma1
###################################################
aa=Box.test.2(residuals(mod2), nlag=c(3,6,9,12),type="Ljung-Box",
decim=4,fitdf=1)
colnames(aa)= c("Retard","p-val.")
t(aa)


###################################################
### code chunk number 33: Anx4.Rnw:581-588
###################################################
theo =  TacvfARMA(theta=c(.3,-.6),lag.max=20)
(acf.theo = theo[-1]/theo[1])
(pacf.theo  = PacfDL(theo/theo[1], LinearPredictor=TRUE)$Pacf)
set.seed(12)
y=arima.sim(n=200,list(ma=c(-0.3, .6)),sd = sqrt(1.5))
(acf.emp= acf(y, 20, plot=FALSE)$acf[-1])
(pacf.emp= pacf(y, 20, plot=FALSE)$acf)


###################################################
### code chunk number 34: Anx4.Rnw:606-614
###################################################
theo =  TacvfARMA(phi=-.8,lag.max=20)
(acf.theo = theo[-1]/theo[1])
(pacf.theo  = PacfDL(theo/theo[1], LinearPredictor=TRUE)$Pacf)
set.seed(23)
y=arima.sim(n=200,list(ar=-.8),
sd = sqrt(1.5))-10
(acf.emp= acf(y, 20, plot=FALSE)$acf[-1])
(pacf.emp= pacf(y, 20, plot=FALSE)$acf)


###################################################
### code chunk number 35: prep_arma1_21a
###################################################
set.seed(4123)
yc = arima.sim(n=200,list(ar=-0.8, ma= c(-0.3, 0.6)), sd = sqrt(1.5))-10


###################################################
### code chunk number 36: prep_arma1_21a0
###################################################
acf.th = ARMAacf(ar=-0.8, ma=c(-0.3, 0.6),lag.max=20, pacf=FALSE)
pacf.th= ARMAacf(ar=-0.8, ma=c(-0.3, 0.6),lag.max=20, pacf=TRUE)


###################################################
### code chunk number 37: plot.arma1 (eval = FALSE)
###################################################
## plotacfthemp(yc, ar=-0.8, ma=c(-0.3, 0.6), lag.max=20)


###################################################
### code chunk number 38: Anx4.Rnw:659-666
###################################################
.PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
pdf(file=paste(file,".pdf",sep=""), width = 7, height = 7, pointsize = 12, bg = "white")
plotacfthemp(yc, ar=-0.8, ma=c(-0.3, 0.6), lag.max=20)
dev.null <- dev.off()
postscript(file=paste(file,".ps",sep=""), width = 7, height = 7, pointsize = 12, bg = "white",horizontal= FALSE,paper="special")
plotacfthemp(yc, ar=-0.8, ma=c(-0.3, 0.6), lag.max=20)
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 39: estim.12
###################################################
(mod12 = Arima(yc, order=c(1,0,2)))


###################################################
### code chunk number 40: prep_arma1_21b
###################################################
require(TSA)
require(polynom)
set.seed(951)
ya = arima.sim(n = 200, list( ma = c(-0.3, 0.6)), sd = sqrt(1.5))
set.seed(7392)
autopol = polynomial(c(1,0.8))*polynomial(c(1,0,0,0,-0.7))
yd=arima.sim(n=200,list(ar=-autopol[-1],ma=c(0,0.6)),sd=sqrt(1.5))
yd=yd+4


###################################################
### code chunk number 41: id.saiso (eval = FALSE)
###################################################
## res=armasubsets(y=yd,nar=10,nma=10,y.name='yd',ar.method='ols')
## plot(res)


###################################################
### code chunk number 42: Anx4.Rnw:720-727
###################################################
.PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
pdf(file=paste(file,".pdf",sep=""), width = 7, height = 7, pointsize = 12, bg = "white")
res=armasubsets(y=yd,nar=10,nma=10,y.name='yd',ar.method='ols')
plot(res)
dev.null <- dev.off()
postscript(file=paste(file,".ps",sep=""), width = 7, height = 7, pointsize = 12, bg = "white",horizontal= FALSE,paper="special")
res=armasubsets(y=yd,nar=10,nma=10,y.name='yd',ar.method='ols')
plot(res)
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


