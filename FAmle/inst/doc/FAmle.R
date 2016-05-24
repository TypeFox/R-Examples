### R code from vignette source 'FAmle.rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: FAmle.rnw:24-25
###################################################
library(FAmle)


###################################################
### code chunk number 2: FAmle.rnw:44-46 (eval = FALSE)
###################################################
## library(FAmle)
## help(distr)


###################################################
### code chunk number 3: FAmle.rnw:51-59
###################################################
x <- -4:4
x
dnorm(x=x,mean=0,sd=1)
Fx <- pnorm(q=x,mean=0,sd=1)
Fx
qnorm(p=Fx,mean=0,sd=1)
set.seed(123)
rnorm(length(x),mean=0,sd=1)


###################################################
### code chunk number 4: FAmle.rnw:62-68
###################################################
distr(x=x,dist='norm',param=c(0,1),type='d')
Fx <- distr(x=x,dist='norm',param=c(0,1),type='p')
Fx
distr(x=Fx,dist='norm',param=c(0,1),type='q')
set.seed(123)
distr(x=length(x),dist='norm',param=c(0,1),type='r')


###################################################
### code chunk number 5: FAmle.rnw:71-73
###################################################
dnorm(x=x,mean=0,sd=1,log=TRUE)
distr(x=x,dist='norm',param=c(0,1),type='d',log=TRUE)


###################################################
### code chunk number 6: FAmle.rnw:77-82
###################################################
p <- .95
mu <- c(0,1,2)
s <- c(.1,.2,.3)
cbind(R=qnorm(p=p,mean=mu,sd=s),FAmle=distr(x=p,dist='norm',
  param=cbind(mu,s),type='q'))


###################################################
### code chunk number 7: FAmle.rnw:85-88
###################################################
p <- c(.1,.5,.9)
cbind(R=qnorm(p=p,mean=mu,sd=s),FAmle=distr(x=p,dist='norm',
  param=cbind(mu,s),type='q'))


###################################################
### code chunk number 8: FAmle.rnw:101-102 (eval = FALSE)
###################################################
## help(mle)


###################################################
### code chunk number 9: FAmle.rnw:106-108 (eval = FALSE)
###################################################
## help(package='FAmle')
## data(package='FAmle')


###################################################
### code chunk number 10: ex1-1
###################################################
data(yarns)
x <- yarns[,1]
hist(x,freq=F,col='grey',border='white',main='')


###################################################
### code chunk number 11: FAmle.rnw:117-120
###################################################
fit.x <- mle(x=x,dist='weibull',start=c(.1,.1))
fit.x
class(fit.x)


###################################################
### code chunk number 12: ex1-2
###################################################
alpha <- .05
plot(x=fit.x,ci=TRUE,alpha=alpha)


###################################################
### code chunk number 13: FAmle.rnw:143-144
###################################################
names(fit.x)


###################################################
### code chunk number 14: FAmle.rnw:149-152
###################################################
x.obs <- 400
distr(x=x.obs,dist=fit.x$dist,param=fit.x$par.hat,type='p',lower.tail=FALSE)
pweibull(q=x.obs,shape=fit.x$par.hat[1],fit.x$par.hat[2],lower.tail=FALSE)


###################################################
### code chunk number 15: FAmle.rnw:155-156
###################################################
distr(x=x.obs,model=fit.x,type='p',lower.tail=FALSE)


###################################################
### code chunk number 16: FAmle.rnw:161-167 (eval = FALSE)
###################################################
## sorted.S.x <- 1 - fit.x[['x.info']][,'Fz']
## sorted.x <- fit.x[['x.info']][,'z']
## empirical.S.x <- 1 - fit.x[['x.info']][,'Emp']
## plot(sorted.x,sorted.S.x,type='l',col='red',
##   xlab=expression(x),ylab=expression(S(x)==1-F(x)))
## points(sorted.x,empirical.S.x,pch=19,cex=.5)


###################################################
### code chunk number 17: FAmle.rnw:170-173 (eval = FALSE)
###################################################
## hist(x,freq=FALSE,main='',col='cyan')
## fun.x <- function(x) distr(x=x,model=fit.x,type='d')
## curve(fun.x,add=T,col='red')


###################################################
### code chunk number 18: FAmle.rnw:175-187
###################################################
sorted.S.x <- 1 - fit.x[['x.info']][,'Fz']
sorted.x <- fit.x[['x.info']][,'z']
empirical.S.x <- 1 - fit.x[['x.info']][,'Emp']
pdf('FAmle-ex1-3.pdf',width=8,height=4)
layout(matrix(1:2,nr=1))
plot(sorted.x,sorted.S.x,type='l',col='red',
  xlab=expression(x),ylab=expression(S(x)==1-F(x)))
points(sorted.x,empirical.S.x,pch=19,cex=.25)
hist(x,freq=FALSE,main='',col='cyan')
fun.x <- function(x) distr(x=x,model=fit.x,type='d')
curve(fun.x,add=T,col='red')
dev.off()


###################################################
### code chunk number 19: FAmle.rnw:197-200
###################################################
dist.vec <- c('weibull','lnorm','gamma')
fit.3 <- list()
for(i in dist.vec) fit.3[[i]] <- mle(x=x,dist=i,start=c(.1,.1))


###################################################
### code chunk number 20: FAmle.rnw:203-204
###################################################
sapply(fit.3,function(h) c(ad=h$ad,aic=h$aic))


###################################################
### code chunk number 21: ex2-1
###################################################
data(station01AJ010)
y <- station01AJ010
layout(matrix(1:2,nr=1))
plot(as.numeric(names(y)),y,type='b',cex=.4,pch=19,
  cex.axis=.75,xlab='Year',ylab='y')
lines(as.numeric(names(y)),lowess(y)[['y']],col='red')
legend('topleft',inset=.01,col='red',lwd=2,'Lowess',bty='n')
acf(y,main='')


###################################################
### code chunk number 22: ex2-2
###################################################
fit.y <- mle(x=y,dist='gamma',start=c(.1,.1))
fit.y
plot(x=fit.y,ci=TRUE,alpha=alpha)


###################################################
### code chunk number 23: FAmle.rnw:249-253
###################################################
names(fit.y)
cov.hat.y <- solve(fit.y$fit$hessian)
se.asy <- sqrt(diag(cov.hat.y))
se.asy


###################################################
### code chunk number 24: FAmle.rnw:257-261
###################################################
boot.y <- boot.mle(model=fit.y,B=2500)
names(boot.y)
se.boot <- apply(boot.y[['par.star']],2,sd)
cbind(asymptotic.se=se.asy,bootstrap.se=as.numeric(se.boot))


###################################################
### code chunk number 25: ex2-3
###################################################
z1 <- seq(-4*se.asy[1],4*se.asy[1],length.out=100)+fit.y[['par.hat']][1]
z2 <- seq(-4*se.asy[2],4*se.asy[2],length.out=100)+fit.y[['par.hat']][2]
fz12 <- outer(z1,z2,
  function(h1,h2) dmvnorm(cbind(h1,h2),fit.y[['par.hat']],fit.y[['cov.hat']])
)
contour(z1,z2,fz12)
title(xlab=colnames(boot.y[['par.star']])[1])
title(ylab=colnames(boot.y[['par.star']])[2])
points(boot.y[['par.star']],cex=.1,col='red',pch=19)


###################################################
### code chunk number 26: FAmle.rnw:294-295
###################################################
boot.y[['p.value']]


###################################################
### code chunk number 27: FAmle.rnw:301-304
###################################################
RP <- c(2,10,20,50,100,500)
p <- 1 - 1/RP
p


###################################################
### code chunk number 28: FAmle.rnw:307-309
###################################################
Q.hat <- distr(x=p,model=fit.y,type='q')
Q.hat


###################################################
### code chunk number 29: FAmle.rnw:320-324
###################################################
Q.conf.int(p=p,model=fit.y,alpha=alpha,ln=FALSE)[c(1,3),]
Q.conf.int(p=p,model=fit.y,alpha=alpha,ln=TRUE)[c(1,3),]
Q.boot.ci(p=p,boot=boot.y,alpha=alpha)[['percentile']][c(1,3),]
Q.boot.ci(p=p,boot=boot.y,alpha=alpha)[['reflexion']][c(1,3),]


###################################################
### code chunk number 30: ex2-4
###################################################
Q.boot.dist <- sapply(as.list(p),function(h) distr(h,dist=fit.y[['dist']],
  param=boot.y[['par.star']],type='q'))
Q.norm <- sapply(as.list(p),function(h) delta.Q(p=h,model=fit.y,ln=FALSE))
layout(matrix(1:6,nc=2))
for(i in 1:ncol(Q.boot.dist)) {
  hist(Q.boot.dist[,i],freq=F,col='gray',border='white',
    main=paste('T = ',RP,sep='')[i],xlab='',cex.axis=.75,las=2)
  curve(distr(x=x,dist='norm',param=Q.norm[,i],type='d'),
    add=TRUE,col='red')
}


###################################################
### code chunk number 31: FAmle.rnw:361-384
###################################################
dGEV <- function(x,shape=1,scale=1,location=0,log=FALSE)
{
	fx <- 1/scale*(1+shape*((x-location)/scale))^(-1/shape-1)*
	  exp(-(1+shape*((x-location)/scale))^(-1/shape))
	if(log) return(log(fx))
	else return(fx)
}
pGEV <- function(q,shape=1,scale=1,location=0,lower.tail=TRUE,log.p=FALSE)
{
	Fx <- exp(-(1+shape*((q-location)/scale))^(-1/shape))
	if(!lower.tail) Fx <- 1 - Fx
	if(log.p) Fx <- log(Fx)
	return(Fx)
}
qGEV <- function(p,shape=1,scale=1,location=0,lower.tail=TRUE,log.p=FALSE)
{
	if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1 - p
	xF <- location+scale/shape*((-log(p))^(-shape)-1)
	return(xF)
}
rGEV <- function(n,shape=1,scale=1,location=0)
	qgev(runif(n),shape,scale,location)


###################################################
### code chunk number 32: FAmle.rnw:389-391
###################################################
fit.y.gev <- mle(x=y,dist='GEV',c(.1,1,1))
fit.y.gev


###################################################
### code chunk number 33: FAmle.rnw:395-397
###################################################
Q.conf.int(p=p,model=fit.y,alpha=alpha)
Q.conf.int(p=p,model=fit.y.gev,alpha=alpha)


###################################################
### code chunk number 34: FAmle.rnw:407-409
###################################################
library(FAdist)
fit.y.ln3 <- mle(x=y,dist='lnorm3',start=c(1,1,1))


###################################################
### code chunk number 35: FAmle.rnw:422-424
###################################################
data(ColesData)
z <- ColesData[,2]


###################################################
### code chunk number 36: ex4-1
###################################################
fit.z <- mle(x=z,dist='gev',start=c(.1,1,1))
plot(x=fit.z,ci=TRUE,alpha=alpha)


###################################################
### code chunk number 37: ex4-2
###################################################
trans.z <- list(function(x) x, function(x) exp(x), function(x) x)
bayes.z <- metropolis(model=fit.z,iter=10000,tun=2,trans.list=trans.z)
plot(x=bayes.z,plot.type='carlin')
bayes.z


###################################################
### code chunk number 38: ex4-3
###################################################
p <- c(.5,.9,.99)
Q.p.post <- sapply(as.list(p),
  function(h) distr(x=h,dist=fit.z[['dist']],
    param=bayes.z[['sims']],type='q'))
layout(matrix(1:length(p),nr=1))
for(i in 1:ncol(Q.p.post)) hist(Q.p.post[,i],
  freq=FALSE,col='steelblue4',
  main=paste('T = ',1/(1-p[i]),sep=''),xlab='')


###################################################
### code chunk number 39: FAmle.rnw:474-478 (eval = FALSE)
###################################################
## prior.z.weak <- function(x)
##   dnorm(x[1],0,1000)*dnorm(x[2],0,1000)*dnorm(x[3],0,1000)
## bayes.z.weak <- metropolis(model=fit.z,iter=10000,
##   tun=2,trans.list=trans.z,prior=prior.z.weak)


###################################################
### code chunk number 40: FAmle.rnw:482-486 (eval = FALSE)
###################################################
## prior.z.weak.2 <- function(x) dmvnorm(x=x[1:3],
##   mean=rep(0,3),sigma=diag(1000,3))
## bayes.z.weak.2 <- metropolis(model=fit.z,iter=10000,
##   tun=2,trans.list=trans.z,prior=prior.z.weak.2)


###################################################
### code chunk number 41: FAmle.rnw:492-496
###################################################
data(floodsNB)
length(floodsNB)
names(floodsNB)[1:5]
floodsNB[[3]]


###################################################
### code chunk number 42: FAmle.rnw:499-506
###################################################
aucoin.2011 <- names(which(sapply(floodsNB,
  function(i) i[['Aucoin.2011']])))
aucoin.2011
length(aucoin.2011)
NB <- list()
for(i in names(floodsNB))
  if(is.element(i,aucoin.2011)) NB[[i]] <- floodsNB[[i]]


###################################################
### code chunk number 43: FAmle.rnw:516-517
###################################################
st <- '01AJ007'


###################################################
### code chunk number 44: FAmle.rnw:521-526
###################################################
st <- '01AJ007'
station.x <- NB[[st]]
ln.drain.x <- station.x[['ln.drain']]
model.x <- list(x=station.x[['peak']],dist='gev')
model.x


###################################################
### code chunk number 45: FAmle.rnw:530-531
###################################################
trans.list.x <- list(function(x) x, function(x) exp(x), function(x) x)


###################################################
### code chunk number 46: FAmle.rnw:535-536
###################################################
station.x[['ln.drain']]


###################################################
### code chunk number 47: FAmle.rnw:541-567
###################################################
gev.inits <- function(x) {
  sig.hat <- sqrt(var(x)*6/pi^2)
  mu.hat <- mean(x) - digamma(1)*sig.hat
  return(c(0.01,sig.hat,mu.hat))
}
nb.stations <- nb.stations.fits <- list()
for(i in names(NB)) {
  if(substr(i,3,3)=='A' & i != st) {
    temp <- NULL
    try(temp <- mle(NB[[i]][['peak']],'gev',gev.inits(NB[[i]][['peak']])),
      silent=TRUE)
    if(!is.null(temp)) {
      nb.stations[[i]] <- NB[[i]]
      nb.stations.fits[[i]] <- temp
    }
  }
}
f <- function(x) c(x[1],log(x[2:3]))
nb.stations.par.hat <- t(sapply(nb.stations.fits,
  function(i) f(i[['par.hat']])))
colnames(nb.stations.par.hat) <- c('shape','ln.scale','ln.location')
nb.stations.par.hat[1:5,]
nb.stations.ln.drain <- sapply(nb.stations,function(i) i[['ln.drain']])
nb.stations.regs <- lapply(list('shape','ln.scale','ln.location'),
  function(i) lm(nb.stations.par.hat[,i]~nb.stations.ln.drain))
names(nb.stations.regs) <- colnames(nb.stations.par.hat)


###################################################
### code chunk number 48: ex4-4
###################################################

par.old <- par(no.readonly=TRUE)
par(oma=c(4,4,4,4),mar=c(0,4,0,2))
layout(matrix(1:3,nr=3))
plot(nb.stations.ln.drain,nb.stations.par.hat[,'shape'],
  type='n',axes=FALSE,xlab='',ylab='')
abline(nb.stations.regs[[1]],col='red')
points(nb.stations.ln.drain,nb.stations.par.hat[,'shape'],pch=19,cex=.8)
title(ylab='shape',cex.lab=1.2);box()
axis(2,cex.axis=1,las=2)
legend('topright',inset=.01,col='red',lty=c(1,0),lwd=c(1,0),
  c(paste('p-value = ',round(summary(nb.stations.regs[[1]])[['coefficients']][2,4],3),
  sep=''),paste('R-squared = ',round(summary(nb.stations.regs[[1]])[['r.squared']],3))),
  bty='n',cex=1.2)

plot(nb.stations.ln.drain,nb.stations.par.hat[,'ln.scale'],type='n',axes=FALSE,xlab='',ylab='')
abline(nb.stations.regs[[2]],col='red')
points(nb.stations.ln.drain,nb.stations.par.hat[,'ln.scale'],pch=19,cex=.8)
title(ylab='ln(scale)',cex.lab=1.2);box()
axis(4,cex.axis=1,las=2)
legend('topleft',inset=.01,col='red',lty=c(1,0),lwd=c(1,0),
  c(paste('p-value = ',round(summary(nb.stations.regs[[2]])[['coefficients']][2,4],3),sep=''),
   paste('R-squared = ',round(summary(nb.stations.regs[[2]])[['r.squared']],3))),bty='n',cex=1.2)


plot(nb.stations.ln.drain,nb.stations.par.hat[,'ln.location'],type='n',axes=FALSE,xlab='',ylab='')
abline(nb.stations.regs[[3]],col='red')
points(nb.stations.ln.drain,nb.stations.par.hat[,'ln.location'],pch=19,cex=.8)
title(ylab='ln(location)',xlab='',cex.lab=1.2);box()
axis(1,cex.axis=1,las=1);axis(2,cex.axis=1,las=2)
legend('topleft',inset=.01,col='red',lty=c(1,0),lwd=c(1,0),c(paste('p-value = ',
  round(summary(nb.stations.regs[[3]])[['coefficients']][2,4],3),sep=''),
  paste('R-squared = ',round(summary(nb.stations.regs[[3]])[['r.squared']],3))),bty='n',cex=1.2)
title(xlab='ln(drainage)',outer=TRUE,line=2.5,cex.lab=1.3)
par(par.old)


###################################################
### code chunk number 49: FAmle.rnw:625-631
###################################################
mu.x <- sapply(nb.stations.regs,
  function(i) as.numeric(coef(i)[1]+coef(i)[2]*ln.drain.x))
mu.x
residuals.x <- sapply(nb.stations.regs,residuals)
sigma.x <- t(residuals.x)%*%residuals.x/nrow(residuals.x)
sigma.x


###################################################
### code chunk number 50: FAmle.rnw:636-638
###################################################
prior.x <- function(p)
  dmvnorm(x=c(p[1:2],log(p[3])),mean=mu.x,sigma=sigma.x)


###################################################
### code chunk number 51: ex4-5
###################################################
inits.1 <- gev.inits(model.x[['x']])
inits.1[2] <- log(inits.1[2])
bayes.x.1 <- metropolis(model=model.x,iter=5000,
  trans.list=trans.list.x,start=inits.1,variance=diag(.1,3),
  prior=prior.x,pass.down.to.C=TRUE)
plot(bayes.x.1)


###################################################
### code chunk number 52: ex4-6
###################################################
start.x <- bayes.x.1[['M']]
variance.x <- bayes.x.1[['V']]
bayes.x.2 <- metropolis(model=model.x,iter=15000,
  trans.list=trans.list.x,start=start.x,variance=variance.x,
  prior=prior.x,pass.down.to.C=TRUE)
plot(bayes.x.2)


###################################################
### code chunk number 53: ex4-7
###################################################
prior.draws <- rmvnorm(100000,mu.x,sigma.x)
prior.draws[,2:3] <- exp(prior.draws[,2:3])
p <- c(.5,.9,.99)
Q.p.post <- sapply(as.list(p),function(h)
  distr(x=h,dist='gev',param=bayes.x.2[['sims']],type='q'))
Q.p.prior <- sapply(as.list(p),function(h)
  distr(x=h,dist='gev',param=prior.draws,type='q'))
layout(matrix(1:length(p),nr=1))
for(i in 1:ncol(Q.p.post)) {
	hist(Q.p.post[,i],freq=FALSE,col='steelblue4',
	  main=paste('T = ',1/(1-p[i]),sep=''),xlab='')
	lines(density(Q.p.prior[,i],bw=3),col='red')
}


