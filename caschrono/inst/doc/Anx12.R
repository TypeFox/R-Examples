### R code from vignette source 'Anx12.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Anx12.Rnw:135-139
###################################################
owidth <- getOption("width") # largeur des sorties
options(width=60, continue="+ ","warn"=-1 )
.PngNo <- 0
nom.fich = "anx12-bitmap-"


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
### code chunk number 12: bfig5 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 7, height = 4, pointsize = 10, bg = "white")


###################################################
### code chunk number 13: bfigps5 (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""),  width = 7, height =4, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 14: zfig2 (eval = FALSE)
###################################################
## dev.null <- dev.off()


###################################################
### code chunk number 15: zfiginclude (eval = FALSE)
###################################################
## cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 16: sim.manu
###################################################
require(caschrono)
set.seed(1618) ; nobs =300
zt = rnorm(nobs+10)
# initialisation
sig2 = rep(0,nobs+10)
epsil = sig2
y = sig2
for ( i in 2:(nobs+10))
{ sig2[i] = 0.1 + 0.9*epsil[i-1]^2
  epsil[i] = sig2[i]^.5 *zt[i]
  y[i] = 5+epsil[i] }
y = y[-(1:10)]


###################################################
### code chunk number 17: Anx12.Rnw:286-291
###################################################
require(fGarch)
spec.1=garchSpec(model=list(mu=5,omega=0.1,alpha=0.9,beta=0),
   rseed=397)
archsim.1 =   garchSim(extended = TRUE,spec.1, n = 300, n.start=10)
head(archsim.1,2)


###################################################
### code chunk number 18: test.hetero.2
###################################################
moyenne=mean(archsim.1[,1])
ret = c(6,12,18)
a1=Box.test.2(archsim.1[,1],ret,type="Box-Pierce",decim=8)
a2=Box.test.2((archsim.1[,1]-moyenne)^2,ret,type="Box-Pierce",decim=8)
a3=Box.test.2(archsim.1[,3]^2,ret,type ="Box-Pierce",decim=8)
a123 = cbind(a1,a2[,2],a3[,2])
colnames(a123)=
   c("Retard","p.val y1","p.val y1 centre carre","p.val z_t")


###################################################
### code chunk number 19: hetero.tex
###################################################
require(xtable)
xtable(a123 , caption="Test de blancheur, p-values"  , label="blanc.y", digits=4)


###################################################
### code chunk number 20: sim.garch
###################################################
spec=garchSpec(model=list(mu=2,omega=0.09,alpha =c(0.15, 0.3),
  beta = 0.4), rseed=9647)
var.margi = 0.09/(1 - 0.15 - 0.3-0.4)
y = garchSim(spec, n = 420, extended = TRUE)
y1 = y[21:420,1]


###################################################
### code chunk number 21: y1.ext1
###################################################
m1 = mean(y1)
(q5 = quantile(abs(y1-m1), probs=c(.975,.98,.985,.99,.995)))
extrem985 = which(abs(y1-m1) > q5[3])
cat('nombre de points : ',length(extrem985),'\n')
y1b = y1
y1b[extrem985]  = q5[3]* sign(y1[extrem985]-m1)


###################################################
### code chunk number 22: y1.ext2
###################################################
mod1=garchFit(~garch(2,1),data=y1b,trace=FALSE,include.mean=TRUE)


###################################################
### code chunk number 23: Anx12.Rnw:394-396
###################################################
param= c( "alpha2","beta1")
estim = as.matrix(coef(mod1)[param])


###################################################
### code chunk number 24: Anx12.Rnw:402-404
###################################################
# matrice des covariances des estimateurs
vcov.par = mod1@fit$cvar[param,param]


###################################################
### code chunk number 25: test.par.2
###################################################
t.alpha =  (coef(mod1)["alpha2"]-0.3)/mod1@fit$se.coef["alpha2"]
p.val = 2* pnorm(-abs(t.alpha))
cat("p-value : ",p.val,"\n")


###################################################
### code chunk number 26: test.par.1
###################################################
# moyenne th\'eorique
theo = as.matrix(c(0.3, 0.4))
# difference
differ = estim-theo
# stat de test
khi = t(differ)%*%solve(vcov.par)%*%differ;
# statistique de test et p-value
cat("stat de test et p-value : ",c(khi, 1-pchisq(khi,df=2)),"\n")


###################################################
### code chunk number 27: lor.1a
###################################################
require(fBasics)
data(csdl)
aa = returns(csdl, percentage = TRUE)
aab = aa[complete.cases(aa) == TRUE,]
r.csdl = its(aab, as.POSIXct(row.names(aab)))
r.lor <- rangeIts(r.csdl[,"L_Oreal"], start= "2007-12-28")
r.lor.0 = r.lor[1:(length(r.lor)-51)]
r.lor.1 = r.lor[(length(r.lor)-50):length(r.lor)]
mod.r.lor=Arima(r.lor.0@.Data,order=c(0,0,4),include.mean=FALSE,
                fixed=c(NA,NA,0,NA))


###################################################
### code chunk number 28: lor.12
###################################################



###################################################
### code chunk number 29: esti.var
###################################################
var.marg.est<-function(mod)
{param.estim = mod@fit$par
 std.estim = mod@fit$se.coef
 k<-which(names(param.estim)=="omega")
 value = param.estim[k]/(1-sum(param.estim[(k+1):length(param.estim)]))
cat("variance marginale : ",value,"\n")
}
mod.garch.lor=garchFit(~garch(1,1),data=r.lor.0@.Data,trace=FALSE,
   include.mean=TRUE,na.action=na.pass)
summary(mod.garch.lor)
var.marg.est(mod.garch.lor)


###################################################
### code chunk number 30: stock.lor
###################################################
npred=50
pred.arima.lor=predict(mod.r.lor,n.ahead=npred)
pred.garch.lor=predict(mod.garch.lor,n.ahead=npred)

str(pred.arima.lor)
demi = qnorm(0.95)*pred.arima.lor$se
binf.arima.lor =  pred.arima.lor$pred-demi
bsup.arima.lor =  pred.arima.lor$pred+demi

demi = qnorm(0.95)*pred.garch.lor$standardDeviation
binf.garch.lor =  pred.garch.lor$meanForecast-demi
bsup.garch.lor =  pred.garch.lor$meanForecast+demi


###################################################
### code chunk number 31: ajust (eval = FALSE)
###################################################
## par(oma=rep(0.5,4))
## mat.lor = cbind(binf.arima.lor,bsup.arima.lor, 
##  binf.garch.lor,bsup.garch.lor,r.lor.1[1:npred])
## matplot(1:npred,mat.lor,type='l', col='black', lty=c(1,1,2,2,3),
##   lwd=2,xlab="horizon 50", ylab="rendement")
## leg.txt = c("GARCH","AR","realisation")
## legend(14,3, leg.txt, lty=c(1,2,3))


###################################################
### code chunk number 32: previreal.Loreal
###################################################
.PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
pdf(file=paste(file,".pdf",sep=""), width = 7, height = 4, pointsize = 10, bg = "white")
par(oma=rep(0.5,4))
mat.lor = cbind(binf.arima.lor,bsup.arima.lor, 
 binf.garch.lor,bsup.garch.lor,r.lor.1[1:npred])
matplot(1:npred,mat.lor,type='l', col='black', lty=c(1,1,2,2,3),
  lwd=2,xlab="horizon 50", ylab="rendement")
leg.txt = c("GARCH","AR","realisation")
legend(14,3, leg.txt, lty=c(1,2,3))
dev.null <- dev.off()
postscript(file=paste(file,".ps",sep=""),  width = 7, height =4, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")
par(oma=rep(0.5,4))
mat.lor = cbind(binf.arima.lor,bsup.arima.lor, 
 binf.garch.lor,bsup.garch.lor,r.lor.1[1:npred])
matplot(1:npred,mat.lor,type='l', col='black', lty=c(1,1,2,2,3),
  lwd=2,xlab="horizon 50", ylab="rendement")
leg.txt = c("GARCH","AR","realisation")
legend(14,3, leg.txt, lty=c(1,2,3))
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 33: rdan.0
###################################################
r.dan = rangeIts(r.csdl[,"Danone"], start= "2007-12-28")
r.dan.0 = r.dan[1:(length(r.dan)-50)]
r.dan.1 = r.dan[(length(r.dan)-49):length(r.dan)]


###################################################
### code chunk number 34: armagarch.dan
###################################################
mod.arma.dan=Arima(r.dan.0@.Data, order=c(2,0,0),include.mean= FALSE)
mod.garch.dan=garchFit(~garch(1,1),data=r.dan.0@.Data,trace=FALSE,
               include.mean=FALSE,na.action=na.pass)


###################################################
### code chunk number 35: stock.dan
###################################################
npred=50
pred.arima.dan=predict(mod.arma.dan,n.ahead=npred)
pred.garch.dan=predict(mod.garch.dan,n.ahead=npred)
str(pred.arima.dan)
demi = qnorm(0.95)*pred.arima.dan$se
binf.arima.dan =  pred.arima.dan$pred-demi
bsup.arima.dan =  pred.arima.dan$pred+demi
#
demi = qnorm(0.95)*pred.garch.dan$standardDeviation
binf.garch.dan =  pred.garch.dan$meanForecast-demi
bsup.garch.dan =  pred.garch.dan$meanForecast+demi


###################################################
### code chunk number 36: pred.dan (eval = FALSE)
###################################################
## par(oma=rep(0.5,4))
## mat.dan = cbind(binf.arima.dan,bsup.arima.dan, 
##  binf.garch.dan,bsup.garch.dan,r.dan.1[1:npred])
## matplot(1:npred,mat.dan,type='l', col='black', lty=c(1,1,2,2,3),
##   lwd=2,xlab="horizon 50", ylab="rendement")
## leg.txt = c("GARCH","AR","realisation")
## legend(14,3, leg.txt, lty=c(1,2,3))


###################################################
### code chunk number 37: pred.dan.plot
###################################################
.PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
pdf(file=paste(file,".pdf",sep=""), width = 7, height = 4, pointsize = 10, bg = "white")
par(oma=rep(0.5,4))
mat.dan = cbind(binf.arima.dan,bsup.arima.dan, 
 binf.garch.dan,bsup.garch.dan,r.dan.1[1:npred])
matplot(1:npred,mat.dan,type='l', col='black', lty=c(1,1,2,2,3),
  lwd=2,xlab="horizon 50", ylab="rendement")
leg.txt = c("GARCH","AR","realisation")
legend(14,3, leg.txt, lty=c(1,2,3))
dev.null <- dev.off()
postscript(file=paste(file,".ps",sep=""),  width = 7, height =4, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")
par(oma=rep(0.5,4))
mat.dan = cbind(binf.arima.dan,bsup.arima.dan, 
 binf.garch.dan,bsup.garch.dan,r.dan.1[1:npred])
matplot(1:npred,mat.dan,type='l', col='black', lty=c(1,1,2,2,3),
  lwd=2,xlab="horizon 50", ylab="rendement")
leg.txt = c("GARCH","AR","realisation")
legend(14,3, leg.txt, lty=c(1,2,3))
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


