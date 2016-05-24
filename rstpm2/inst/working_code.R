## package.skeleton(name="rstpm2", path="c:/usr/src/R", force=TRUE, namespace=TRUE, code_files="pm2-3.R")
## Local Windows setup:
## Rtools.bat
## R CMD INSTALL --html "c:/usr/src/R/rstpm2/pkg"
## R CMD build "c:/usr/src/R/rstpm2/pkg"
## R CMD INSTALL --build "c:/usr/src/R/rstpm2/pkg"
## R CMD CHECK "c:/usr/src/R/rstpm2/pkg"
##
## Local Ubuntu setup:
## R CMD INSTALL --html ~/src/R/rstpm2/pkg --library=~/R/x86_64-pc-linux-gnu-library/2.12
## R CMD build ~/src/R/rstpm2/pkg
## R CMD build --binary ~/src/R/rstpm2/pkg
##
## testPackage <- TRUE
## if (testPackage) {
##   require(splines)
##   require(survival)
##   require(bbmle)
## }


require(abind)
X <- matrix(seq(0,1,length=5*10),nrow=10)
beta <- seq(0,1,length=5)
H <- exp(as.vector(X %*% beta))
dHdbeta <- X * H # row=indiv, col=beta

d2Hdbeta2 <- aperm(abind(lapply(1:ncol(X), function(k) X[,k] * X * H),along=3),c(2,3,1))
abind(lapply(1:nrow(X), function(i) (X[i,] %*% t(X[i,])) * H[i]),along=3) -
    aperm(abind(lapply(1:ncol(X), function(k) X[,k] * X * H),along=3),c(2,3,1))

numder <- function(f,x,eps=1e-8) (f(x+eps)-f(x-eps))/2/eps
expit <- function(x) 1/(1+exp(-x))
numder(expit,2)
expit(2)*expit(-2)
numder(dnorm,2)
-dnorm(2)*2

refresh
require(rstpm2)
data(brcancer)

## Stata estimated coef for hormon
## PH:     -.3614357
## PO:     -.474102
## Probit: -.2823338
system.time(print( stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)))
system.time(print(pfit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)))
##
system.time(print( stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,type="PO")))
system.time(print(pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,type="PO")))
##
system.time(print( stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,type="probit")))
system.time(print(pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,type="probit"))) # slow

if (FALSE) {
    debug(pstpm2)
    pfit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,sp=1)
    ## towards the end of the pstpm2 function...
    sum(diag(solve(optimHess(coef(mle2),negllsp,sp=1)) %*% optimHess(coef(mle2),negll0sp,sp=1)))
    sum(diag(solve(optimHess(coef(mle2),negllsp,sp=fit$sp)) %*% optimHess(coef(mle2),negll0sp,sp=fit$sp)))
    negllsp(coef(mle2),sp=1)
    negll0sp(coef(mle2),sp=1)
}

## delayed entry
## Stata estimated coef for hormon (PH): -1.162504
data(brcancer)
brcancer2 <- transform(brcancer,startTime=ifelse(hormon==0,rectime*0.5,0))
## debug(stpm2)
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
      logH.formula=~nsx(log(rectime),df=3,stata=TRUE)))
head(predict(fit,se.fit=TRUE))
## delayed entry and tvc
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
                     logH.formula=~nsx(rectime,df=3),
                     tvc.formula=~hormon:nsx(rectime,df=3,stata=TRUE)))
head(predict(fit,se.fit=TRUE)) 
pstpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2)

## additive model
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
                     logH.formula=~nsx(rectime,df=3),
                     tvc.formula=~hormon:nsx(rectime,df=3,stata=TRUE)))

require(rstpm2)
require(frailtypack)
data(dataAdditive)
system.time(mod2n <- pstpm2(Surv(t1,t2,event)~var1,
                           data=dataAdditive,
                           RandDist="LogN",
                           smooth.formula=~s(t2),
                           cluster=dataAdditive$group, nodes=20))

system.time(mod2nb <- stpm2(Surv(t1,t2,event)~var1,
                           data=dataAdditive,
                           RandDist="LogN",
                           logH.formula=~ns(t2,df=7),
                           cluster=dataAdditive$group, nodes=20))

mod1 <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+var1,data=dataAdditive,
                     n.knots=8,kappa1=0.1,cross.validation=TRUE)
mod1n <- frailtyPenal(Surv(t1,t2,event)~cluster(group)+var1,data=dataAdditive,
                     n.knots=8,kappa1=0.1,cross.validation=TRUE, RandDist="LogN")

system.time(mod2 <- stpm2(Surv(t1,t2,event)~var1, # Gamma
                          data=dataAdditive,
                          logH.formula=~ns(t2,df=7),
                          cluster=dataAdditive$group))

system.time(coxph1 <- coxph(Surv(t1,t2,event)~var1+frailty(group,distribution="gaussian"),
                          data=dataAdditive))
summary(coxph1)

system.time(mod2n <- stpm2(Surv(t1,t2,event)~var1,
                           data=dataAdditive,
                           RandDist="LogN",
                           optimiser="NelderMead",
                           logH.formula=~ns(t2,df=7),
                           cluster=dataAdditive$group, nodes=20))
system.time(mod2nb <- stpm2(Surv(t1,t2,event)~var1,
                           data=dataAdditive,
                           RandDist="LogN",
                           logH.formula=~ns(t2,df=7),
                           cluster=dataAdditive$group, nodes=20))

system.time(mod3 <- pstpm2(Surv(t1,t2,event)~var1,
                           data=dataAdditive,
                           RandDist="LogN",
                           criterion="BIC",
                           smooth.formula=~s(t2),
                           cluster=dataAdditive$group, nodes=20))

system.time(mod3 <- coxph(Surv(t1,t2,event)~frailty(group,distribution="gamma")+var1,data=dataAdditive))
summary(mod2)
coef2 <- coef(summary(mod2))
theta <- exp(coef2[nrow(coef2),1])
se.logtheta <- coef2[nrow(coef2),2]
se.theta <- theta*se.logtheta
test.statistic <- 1/se.logtheta
pchisq(test.statistic,df=1,lower.tail=FALSE)/2

refresh
require(rstpm2)
require(ICE)
data(ICHemophiliac)
ICHemophiliac2 <- transform(as.data.frame(ICHemophiliac),event=3)
fit1 <- stpm2(Surv(left,right,event,type="interval")~1,data=ICHemophiliac2)
estimate <- ickde(ICHemophiliac, m=200, h=0.9)
plot(estimate, type="l", ylim=c(0,0.20))
tt <- seq(0,20,length=301)[-1]
## plot(fit1,newdata=data.frame(x=1),type="density",add=TRUE,line.col="blue")
lines(tt,predict(fit1,newdata=data.frame(right=tt),type="density"),col="blue")

## reg1 <- survreg(Surv(left,right,event,type="interval")~1,data=ICHemophiliac2)
## weibullShape <- 1/reg1$scale
## ## weibullScale <- exp(predict(reg1,type="lp"))
## weibullScale <- predict(reg1);
## tt <- seq(0,20,length=301)
## estimate <- ickde(ICHemophiliac, m=200, h=0.9)
## plot(estimate, type="l", ylim=c(0,0.15))
## lines(tt,dweibull(tt,weibullShape,weibullScale),lty=2)


## two-dimensional smoothers
x1 <- x2 <- seq(0,1,length=11)
dat <- expand.grid(x1=x1,x2=x2)
dat$y <- rnorm(nrow(dat))
require(mgcv)
fit <- gam(y~s(x1,x2),data=dat)
fit$smooth


system.time(print(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,type="probit")))
system.time(print(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,type="probit",use.rcpp=FALSE)))

system.time(print(fit2 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer, type="PO")))
system.time(print(fit2 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer, type="PO",use.rcpp=FALSE)))

system.time(print(stpm2Gen(Surv(rectime,censrec==1)~hormon,data=brcancer)))
system.time(print(stpm2Gen(Surv(rectime,censrec==1)~hormon,data=brcancer, use.rcpp=FALSE)))
head(predict(fit,se.fit=TRUE))
head(predict(fit,type="haz",se.fit=TRUE))

brcancer <- brcancer[rep(1:nrow(brcancer),each=500),]
system.time(print(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer))) # faster than Stata!
system.time(print(pfit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)))
plot(pfit,newdata=data.frame(hormon=0))

refresh
require(rstpm2)
data(brcancer)
system.time(print(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                               tvc=list(hormon=3))))
system.time(print(pfit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer,sp.init=c(0.0001,0.0001),
                                 tvc.formula=~s(log(rectime),by=hormon))))
print(pstpm2(Surv(rectime,censrec==1)~1,data=brcancer,init=coef(pfit)*100,
                                 tvc.formula=~s(log(rectime),by=hormon)))

summary(pfit)
plot(pfit,newdata=data.frame(hormon=0))
plot(pfit,newdata=data.frame(hormon=1),add=TRUE)
plot(pfit,newdata=data.frame(hormon=0),type="haz")
plot(pfit,newdata=data.frame(hormon=1),type="haz",add=TRUE)

pfit <- pstpm2(Surv(rectime/365,censrec==1)~1,data=brcancer) # OK
plot(pfit,newdata=data.frame(hormon=0))
system.time(print(pfit <- pstpm2(Surv(rectime/365,censrec==1)~1,data=brcancer,
               tvc.formula=~s(log(rectime/365),by=hormon))))
plot(pfit,newdata=data.frame(hormon=0)) # OK


times <- seq(500,2000,by=500)
meansurv1 <- t(sapply(times,function(time) predict(pfit,transform(brcancer,rectime=time,hormon=1),type="meansurv",se.fit=TRUE)))
meansurv0 <- t(sapply(times,function(time) predict(pfit,transform(brcancer,rectime=time,hormon=0),type="meansurv",se.fit=TRUE)))
matplot(times,meansurv1,type="l",lty=c(1,2,2),col=1)
matlines(times,meansurv0,type="l",lty=c(1,2,2),col=2)

meansurvdiff1 <- t(sapply(times,function(time) predict(pfit,transform(brcancer,rectime=time,hormon=0),type="meansurvdiff",var="hormon",se.fit=TRUE)))
matplot(times,meansurvdiff1,type="l",lty=c(1,2,2),col=1)



system.time(print(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,control=list(parscale=100,reltol=1e-10),use.rcpp=FALSE)))

summary(fit)
system.time(print(fit2 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,control=list(parscale=10000.0),reltol=1e-10,init=0.0001*coef(fit))))
summary(fit2)
plot(fit2,newdata=data.frame(hormon=1))

brcancerN <- brcancer[rep(1:nrow(brcancer),each=100),]
system.time(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancerN,use.rcpp=FALSE,
                         control=list(parscale=0.1,reltol=1e-10)))
summary(fit)
system.time(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancerN,use.rcpp=TRUE))
summary(fit)

###### penalised likelihood
## environment(pstpm2) <- environment(rstpm2::pstpm2)
## require(rstpm2)
try(detach("package:rstpm2",unload=TRUE))
## source("/home/MEB/marcle/src/R/rstpm2/R/pm2-3.R")

refresh
require(rstpm2)
data(brcancer)
brcancer$recyear <- brcancer$rectime/365
system.time(fit0 <- stpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,df=5))
system.time(pfit0 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,sp.init=1))
system.time(pfit0.1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                smooth.formula=~s(log(recyear),k=15),sp.init=10,alpha=2))
system.time(pfit1.1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                smooth.formula=~s(log(recyear)),sp.init=10,criterion="BIC"))
system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                smooth.formula=~s(recyear),sp.init=10))


plot(pfit0,newdata=data.frame(hormon=1),line.col="red",type="hazard")
plot(pfit0.1,newdata=data.frame(hormon=1),line.col="blue",add=TRUE,type="hazard")
plot(pfit1.1,newdata=data.frame(hormon=1),line.col="orange",add=TRUE,type="hazard")
plot(fit0,newdata=data.frame(hormon=1),line.col="green",type="hazard",add=TRUE)
plot(pfit2,newdata=data.frame(hormon=1),line.col="black",type="hazard",add=TRUE)

plot(pfit0,newdata=data.frame(hormon=1),line.col="red")
plot(pfit0.1,newdata=data.frame(hormon=1),line.col="blue",add=TRUE)
plot(pfit1.1,newdata=data.frame(hormon=1),line.col="pink",add=TRUE)
plot(fit0,newdata=data.frame(hormon=1),line.col="green",add=TRUE)
plot(pfit2,newdata=data.frame(hormon=1),line.col="black",add=TRUE)



system.time(pfit0.check <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer, sp=pfit0@sp, use.rcpp=FALSE))
system.time(pfit0.check2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer, sp=pfit0@sp))
summary(pfit0)
summary(pfit0.check)
summary(pfit0.check2)

system.time(pfit0 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=30),sp.init=1))
system.time(pfit0.1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=30),sp.init=1,criterion="BIC"))
plot(pfit0,newdata=data.frame(hormon=1),line.col="red",type="hazard")
plot(pfit0.1,newdata=data.frame(hormon=1),line.col="blue",add=TRUE,type="hazard")

system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=30),sp=10,pen="h",
                            smoother.parameters=list("log(recyear)"=list(var="recyear",
                                                         inverse=exp,
                                                         transform=log))))
plot(pfit1,newdata=data.frame(hormon=1),line.col="green",add=TRUE,type="hazard")


system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=30),sp=10,criterion="GCV"))
system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=20),sp=10,criterion="BIC"))
plot(pfit1,newdata=data.frame(hormon=1),type="hazard",ylim=c(0,0.25))
plot(pfit2,newdata=data.frame(hormon=1),add=TRUE,line.col="blue",type="hazard")

system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=30),sp=1,criterion="GCV"))

system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(recyear,k=20),sp=1,use.rcpp=F,penalty="h"))
system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(recyear,k=20),sp=1,penalty="h"))
plot(pfit2,newdata=data.frame(hormon=1),type="hazard")

system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(recyear,k=30),sp=1,use.rcp=FALSE))

plot(pfit1,newdata=data.frame(hormon=1),type="hazard")
plot(pfit2,newdata=data.frame(hormon=1),type="hazard")

system.time(pfit2.0 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(recyear,k=30),sp=0.055,penalty="h",cr="GCV"))
system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(recyear,k=30),sp=0.055,use.rcp=FALSE,penalty="h"))
rstpm2:::gcv(pfit2)
plot(pfit2,newdata=data.frame(hormon=1),line.col="red",add=TRUE,type="hazard")
plot(pfit2.0,newdata=data.frame(hormon=1),line.col="green",add=TRUE,type="hazard")

require(frailtypack)
fpack1 <- frailtyPenal(Surv(recyear,censrec==1)~hormon, data=brcancer, cross.validation=TRUE, n.knots=10, kappa1=0.1)
plot(fpack1)

system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon+x3,data=brcancer,
                logH.formula=~s(log(recyear),k=20)+s(x3),sp=c(0.1,0.1)))
system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon+x3,data=brcancer,
                logH.formula=~s(log(recyear),k=20)+s(x3)))
system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=20)+s(x3)))
plot(pfit1,newdata=data.frame(hormon=1,x3=20))
plot(pfit1,newdata=data.frame(hormon=0,x3=20),type="hazard")
plot(pfit1,newdata=data.frame(hormon=1,x3=20),type="hazard",add=TRUE,line.col="blue",lty=1)

summary(pfit1)
brcancerN <- brcancer[rep(1:nrow(brcancer),each=100),]
system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancerN,
                logH.formula=~s(log(recyear),k=20)))
plot(pfit1,newdata=data.frame(hormon=1))
pfit1@gam$sp
par(mfrow=c(2,2))
plot(pfit1,newdata=data.frame(hormon=1))
summary(pfit1@gam)$edf
rstpm2:::gcv(pfit1)
rstpm2:::gcvc(pfit1,nn)
sps <- as.list(10^(seq(-4,2,by=0.5)))
system.time(pfit2 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                logH.formula=~s(log(recyear),k=20), sp=sps))
gcvs <- lapply(pfit2,rstpm2:::gcv)
plot(sps,unlist(gcvs),type="l",log="x")
plot(sapply(gcvs,attr,"negll"),sapply(gcvs,attr,"trace"),type="l",asp=1)
plot(sapply(gcvs,attr,"trace"),sapply(gcvs,attr,"negll"),type="l",asp=1)
plot(sps,sapply(pfit2,rstpm2:::aicc,nn=nn),type="l",log="x")
plot(sps,sapply(pfit2,rstpm2:::bic,nn=nn),type="l",log="x")
##gcvc
brcancer$recyear <- brcancer$rectime/365
sps <- 10^(seq(-4,2,by=0.5))
gcvcs <- sapply(sps, function(sp) {
gcvc(pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                      logH.formula=~s(recyear,k=30), sp=sp),length(brcancer$recyear))
})
plot(sps,gcvcs,type="l",log="x")
###bic
brcancer$recyear <- brcancer$rectime/365
sps <- 10^(seq(-4,2,by=0.5))
gcvcs <- sapply(sps, function(sp) {
  bic(pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
              logH.formula=~s(recyear,k=30), sp=sp),length(brcancer$recyear))
})
plot(sps,gcvcs,type="l",log="x")
###aicc
brcancer$recyear <- brcancer$rectime/365
sps <- 10^(seq(-4,2,by=0.5))
gcvcs <- sapply(sps, function(sp) {
  aicc(pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
             logH.formula=~s(recyear,k=30), sp=sp),length(brcancer$recyear))
})
plot(sps,gcvcs,type="l",log="x")
#########################


### penalty functions
require(mgcv)
require(gaussquad)

## Outline:
## get w, lambda, X0, X1, X2, X3
## calculate s0, s1, s2, s3
## calculate h2 and pfun=integrate(h2^2,t)
## calculate dh2sq.dbeta and dpfun=integrate(dh2sq.dbeta,t)
##
## calculate w, lambda, X0, X1, X2, X3
derivativeDesign <- 
function (functn, lower = -1, upper = 1, rule = NULL,
    ...) 
{
    pred <- if (length(list(...)) && length(formals(functn)) > 
              1) 
        function(x) functn(x, ...)
    else functn
    if (is.null(rule))
        rule <-    ## gaussquad::legendre.quadrature.rules(20)[[20]]
        data.frame(x = c(0.993128599185095, 0.963971927277914, 0.912234428251326, 
                       0.839116971822219, 0.746331906460151, 0.636053680726515, 0.510867001950827, 
                       0.37370608871542, 0.227785851141646, 0.0765265211334977, -0.0765265211334974, 
                       -0.227785851141645, -0.373706088715418, -0.510867001950827, -0.636053680726516, 
                       -0.746331906460151, -0.839116971822219, -0.912234428251326, -0.963971927277913, 
                       -0.993128599185094),
                   w = c(0.0176140071391522, 0.040601429800387, 
                       0.0626720483341092, 0.0832767415767053, 0.101930119817241, 0.11819453196152, 
                       0.131688638449176, 0.14209610931838, 0.149172986472603, 0.152753387130726, 
                       0.152753387130726, 0.149172986472603, 0.142096109318381, 0.131688638449175, 
                       0.11819453196152, 0.10193011981724, 0.0832767415767068, 0.0626720483341075, 
                       0.0406014298003876, 0.0176140071391522))
    lambda <- (upper - lower)/(2)
    mu <- (lower + upper)/(2)
    x <- lambda * rule$x + mu
    w <- rule$w
    eps <- .Machine$double.eps^(1/8)
    X0 <- pred(x)
    X1 <- (-pred(x+2*eps)+8*pred(x+eps)-8*pred(x-eps)+pred(x-2*eps))/12/eps
    X2 <- (-pred(x+2*eps)/12+4/3*pred(x+eps)-5/2*pred(x)+4/3*pred(x-eps)-pred(x-2*eps)/12)/eps/eps
    X3 <- (-pred(x+3*eps)/8+pred(x+2*eps)-13/8*pred(x+eps)+
           13/8*pred(x-eps)-pred(x-2*eps)+pred(x-3*eps)/8)/eps/eps/eps
    return(list(x=x,w=w,lambda=lambda,X0=X0,X1=X1,X2=X2,X3=X3))
}
hpfun <- function(beta,design) {
    lapply(design,function(obj) {
        s0 <- as.vector(obj$X0 %*% beta)
        s1 <- as.vector(obj$X1 %*% beta)
        s2 <- as.vector(obj$X2 %*% beta)
        s3 <- as.vector(obj$X3 %*% beta)
        h2 <- (s3+3*s1*s2+s1^3)*exp(s0)
        obj$lambda*sum(obj$w*h2^2)
    })
}
hdpfun <- function(beta,design) {
    lapply(design, function(obj) {
        s0 <- as.vector(obj$X0 %*% beta)
        s1 <- as.vector(obj$X1 %*% beta)
        s2 <- as.vector(obj$X2 %*% beta)
        s3 <- as.vector(obj$X3 %*% beta)
        h2 <- (s3+3*s1*s2+s1^3)*exp(s0)
        dh2sq.dbeta <- 2*h2*(exp(s0)*(obj$X3+3*(obj$X1*s2+obj$X2*s1)+3*s1^2*obj$X1)+h2*obj$X0)
        obj$lambda*colSums(obj$w*dh2sq.dbeta)
    })
}
smootherDesign <- function(gamobj,data) {
    d <- data[1,,drop=FALSE] ## how to get mean prediction values, particularly for factors?
    makepred <- function(var) {
        function(value) {
            d <- d[rep(1,length(value)),]
            d[[var]] <- value
            predict(gamobj,newdata=d,type="lpmatrix")
        }
    }
    lapply(gamobj$smooth, function(smoother) {
        var <- smoother$term
        pred <- makepred(var)
        derivativeDesign(pred,lower=min(data[[var]]),upper=max(data[[var]]))
    })
}
## example data
d <- within(data.frame(x=seq(0,1,length=301)), {
    mu <- exp(x)
    y <- rnorm(301,mu,0.01)
})
fit <- gam(y~s(x),data=d,family=gaussian(link="log"))
beta <- coef(fit)
design <- smootherDesign(fit,d)
hpfun(beta,design)
hdpfun(beta,design)

    

## Testing...
require(mgcv)
d <- within(data.frame(x=seq(0,2*pi,length=301)), {
    mu <- sin(x)
    dmu <- cos(x)
    y <- rnorm(301,mu,0.001)
})
fit <- gam(y~s(x),data=d)
mat <- predict(fit,type="lpmatrix")
with(d,plot(x,y))
with(d,lines(x,mu,lwd=2))
with(d,lines(x,predict(fit),col="blue",lwd=2))
par(mfrow=c(3,2))
pred <- function(eps,obj=fit,data=d,var="x") {
    nd <- d
    nd[[var]] <- nd[[var]]+eps
    predict(obj,newdata=nd,type="lpmatrix")
}
## First derivative
eps <- .Machine$double.eps^(1/8)
matD <- (pred(eps) - pred(-eps)) / 2 / eps
with(d,plot(x,dmu,lwd=2,type="l"))
with(d,lines(x,matD %*% coef(fit),col="blue",lwd=2))
##
## 1/12 	−2/3 	0 	2/3 	−1/12
eps <- .Machine$double.eps^(1/8)
matD <- (-pred(2*eps)+8*pred(eps)-8*pred(-eps)+pred(-2*eps))/12/eps
with(d,plot(x,dmu,lwd=2,type="l"))
with(d,lines(x,matD %*% coef(fit),col="blue",lwd=2))
##
## Second derivative
eps <- .Machine$double.eps^(1/8)
matD2 <- (pred(eps)-2*pred(0)+pred(-eps))/eps/eps
with(d,plot(x,-mu,lwd=2,type="l"))
with(d,lines(x,matD2 %*% coef(fit),col="blue",lwd=2))
##
## −1/12 	4/3 	−5/2 	4/3 	−1/12 	
eps <- .Machine$double.eps^(1/8)
matD2 <- (-pred(2*eps)/12+4/3*pred(eps)-5/2*pred(0)+4/3*pred(-eps)-pred(-2*eps)/12)/eps/eps
with(d,plot(x,-mu,lwd=2,type="l"))
with(d,lines(x,matD2 %*% coef(fit),col="blue",lwd=2))
##
## Third derivatives
eps <- .Machine$double.eps^(1/8)
matD3 <- (pred(2*eps)-
          2*pred(eps)+
          2*pred(-eps)-
          pred(-2*eps))/2/eps/eps/eps
with(d,plot(x,-dmu,lwd=2,type="l"))
with(d,lines(x,matD3 %*% coef(fit),col="blue",lwd=2))
##
## 1/8 	−1 	13/8 	0 	−13/8 	1 	−1/8
eps <- .Machine$double.eps^(1/8)
matD3 <- (-pred(3*eps)/8+pred(2*eps)-13/8*pred(eps)+
          13/8*pred(-eps)-pred(-2*eps)+pred(-3*eps)/8)/eps/eps/eps
with(d,plot(x,-dmu,lwd=2,type="l"))
with(d,lines(x,matD3 %*% coef(fit),col="blue",lwd=2))


## (browse-url "http://en.wikipedia.org/wiki/Finite_difference_coefficients")



require(mgcv)
data <- data.frame(x=1:10,y=1:10)
fit <- gam(y~s(x,k=5,bs="ps"),data=data)

round(cbind(1,(spline.des(knots=fit$smooth[[1]]$knots,x=data$x)$design %*% 
  qr.Q(attr(fit$smooth[[1]],"qrc"),complete=TRUE))[,-1]) -
           predict(fit,type="lpmatrix"),1e-10)

round(cbind(1,(spline.des(knots=fit$smooth[[1]]$knots,x=5:6)$design %*% 
                 qr.Q(attr(fit$smooth[[1]],"qrc"),complete=TRUE))[,-1]) -
        predict(fit,newdata=data.frame(x=5:6),type="lpmatrix"),1e-10)

cbind(0,(spline.des(knots=fit$smooth[[1]]$knots,x=data$x,derivs=rep(1,nrow(data)))$design %*% 
        qr.Q(attr(fit$smooth[[1]],"qrc"),complete=TRUE))[,-1]) -
  (predict(fit,newdata=transform(data,x=x+1e-5),type="lpmatrix")-
     (predict(fit,newdata=transform(data,x=x-1e-5),type="lpmatrix")))/2e-5



#######Optimal fitting#######
###GCV,AICC,BIC or GCVC to choose smoothing parameters###
opt.val<-function(pstpm2.fit,nn){
  like<-pstpm2.fit@like
  Hl<-numDeriv::hessian(like,coef(pstpm2.fit))
  Hinv<-vcov(pstpm2.fit)
  trace<-sum(diag(Hinv%*%Hl))
  loglike<-(like(coef(pstpm2.fit)))/nn
  gcv<-(trace-loglike)/nn
  aicc<-(-2*loglike+2*trace*nn/(nn-trace-1))/nn
  bic<-(-2*loglike+trace*log(nn))/nn
  gcvc<-(-2*loglike-2*nn*log(1-trace/nn))/nn
  out<-c(loglike,gcv,aicc,bic,gcvc)
  return(out)
}
###############################
###############################
# setClass("opt.fit", representation(
#                                    num.ind = "numeric",
#                                    cr = "numeric",
#                                    tops = "data.frame",
#                                    sp.opt = "numeric",
#                                    fun.min = "numeric"
# ),
#          contains="pstpm2")
# #########################
# opt.fit<-function(formula,data,logH.formula,sp.low,sp.upp,num.sp,timeVar = NULL){
#   ###number of individual
#   num.ind <- nrow(data)
#   #####Censoring rate####
#   ## set up the data
#   ## ensure that data is a data frame
#   data <- get_all_vars(formula, data)
#   #   ## parse the function call
#   #   Call <- match.call()
#   #   mf <- match.call(expand.dots = FALSE)
#   #   m <- match(c("formula", "data", "subset", "contrasts", "weights"),
#   #              names(mf), 0L)
#   #   mf <- mf[c(1L, m)]
#   stopifnot(length(lhs(formula))>=2)
#   eventExpr <- lhs(formula)[[length(lhs(formula))]]
#   delayed <- length(lhs(formula))==4
#   timeExpr <- lhs(formula)[[if (delayed) 3 else 2]] # expression
#   if (is.null(timeVar))
#     timeVar <- all.vars(timeExpr)
#   time <- eval(timeExpr, data)
#   if (delayed) {
#     time0Expr <- lhs(formula)[[2]]
#     time0 <- eval(time0Expr, data)
#   }
#   event <- eval(eventExpr,data)
#   cr <- sum(event > min(event))/num.ind
#   #   
#   #   cr=table(lhs(formula)[[if (delayed) 4 else 3]][2])/nn
#   ##nn<-length(brcancer$recyear)
#   # system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
#   #                             logH.formula=~s(recyear,k=30), sp=1e-1))
#   # plot(pfit1,newdata=data.frame(hormon=1))
#   
#   #sps <- 10^(seq(-4,4,by=0.5))
#   #   sp.low=10^-4
#   #   sp.upp=4000
#   #   num.sp=30
#   sps <- 10^(seq(log10(sp.low),log10(sp.upp),length=num.sp))
#   optvals <- sapply(sps, function(sp) {
#     opt.val(pstpm2(formula,data,logH.formula=NULL, sp=sp),num.ind)
#   })
#   tops<-t(optvals)
#   colnames(tops) <- c("loglike","gcv","aicc","bic","gcvc")
#   rownames(tops) <- rownames(tops, do.NULL = FALSE, prefix = "Obs.")
#   # tops<-as.data.frame(tops)
#   tops<-as.data.frame(tops)
#   ####Plot#########
#   #par(mfrow=c(1,2))
#   ###to choose optimal smoothing parameter ###
#   ind.min <- sapply(2:5,function(x) order(tops[,x])[1])
#   sp.opt <- sps[ind.min]
#   obj<-pstpm2(formula,data,logH.formula=NULL, sp=sp.opt[1])
#   fun.min <- sapply(2:5,function(x) min(tops[,x]))
#   # if(ind.min[1]==1)
#   # stop("Hit left boundary, make sp.low smaller.")
#   # if(ind.min[1]==num.sp)
#   # stop("Hit right boundary, make sp.upp bigger.")
#   #   with(tops,matplot(sps,tops[,-1],type="l",col=1:4,lty=1:4,xlab="x",ylab="y"))
#   #   points(sp.opt,fun.min,pch=4,lwd=2,cex=1.2)
#   #   lines(sp.opt,fun.min,err=-1,col=1:4,lty=1:4)
#   
#   ###Estimate final model with optimal value of sp###
#   
#   #   
#   #   summary(pfit.obj)
#   #########################################
#   out <- as(obj,"opt.fit")
#   out <- new("opt.fit",
#              coef = pstpm2@coef,
#              fullcoef = pstpm2@fullcoef,
#              vcov = pstpm2@vcov,
#              min = pstpm2@min,
#              details = pstpm2@details,
#              minuslogl = pstpm2@minuslogl,
#              method = pstpm2@method,
#              data = data,
#              formula = pstpm2@formula,
#              optimizer = "optim",
#              xlevels = .getXlevels(pstpm2@terms,pstpm2@model.frame),
#              ##contrasts = attr(X, "contrasts"),
#              contrasts = NULL, # wrong!
#              logli = pstpm2@logli,
#              ##weights = weights,
#              Call = pstpm2@Call,
#              terms = pstpm2@terms,
#              model.frame = pstpm2@model.frame,
#              gam = pstpm2@gam,
#              timeVar = pstpm2@timeVar,
#              timeExpr = pstpm2@timeExpr,
#              like = pstpm2@like,
#              negll<-pstpm2@negll,
#              call.formula = pstpm2@call.formula,
#              x = pstpm2@x,
#              xd = pstpm2@xd,
#              termsd = pstpm2@termsd, # wrong!
#              y = pstpm2@y,
#              num.ind = num.ind,
#              cr = cr,
#              tops = tops,
#              sp.opt = sp.opt,
#              fun.min = fun.min)
#   
#   return(out)
# }


#####load data####
load("brcancer.rda")
data(brcancer)
brcancer$recyear <- brcancer$rectime/365
####model fit###
opt.fit(Surv(recyear,censrec==1)~hormon,data=brcancer,
        logH.formula=~s(recyear), sp.low=10^-4,sp.upp=4000,
        num.sp=30,timeVar = NULL)
# ###methods for Plot ###
# setMethod(
#   f= "plot",
#   signature(x="opt.fit", y="missing"),
#   definition=function (x,y,...){
#     matplot(x@sps,x@tops[,-1],type="l",col=1:4,lty=1:4,xlab="",ylab="")
#     points(x@sp.opt,x@fun.min,pch=4,lwd=2,cex=1.2)
#     lines(x@sp.opt,x@fun.min,err=-1,col=1:4,lty=1:4)
#   }
# )
# ####methods for print####
# setMethod ("print",signature(x="opt.fit", y="missing"),
#            function(x,...){
#              cat("*** Class opt.fit, method Print *** \n")
#              cat("* Optimal SP ="); print (x@sp.opt)
#              cat("* GCV = \n"); print (x@fun.min[1])
#              cat("******* End Print (opt.fit) ******* \n")
#            }
# )
##########################

aplot <- 
function (x, y, ...) 
{
    .local <- function (x, y, newdata, type = "surv", xlab = NULL, 
        ylab = NULL, line.col = 1, ci.col = "grey", lty = par("lty"), 
        add = FALSE, ci = !add, rug = !add, var = NULL, ...) 
    {
        browser()
        y <- predict(x, newdata, type = type, var = var, grid = TRUE, 
            se.fit = TRUE)
        if (is.null(xlab)) 
            xlab <- deparse(x@timeExpr)
        if (is.null(ylab)) 
            ylab <- switch(type, hr = "Hazard ratio", hazard = "Hazard", 
                surv = "Survival", sdiff = "Survival difference", 
                hdiff = "Hazard difference", cumhaz = "Cumulative hazard")
        xx <- attr(y, "newdata")
        xx <- eval(x@timeExpr, xx)
        if (!add) 
            matplot(xx, y, type = "n", xlab = xlab, ylab = ylab, 
                ...)
        if (ci) 
            polygon(c(xx, rev(xx)), c(y[, 2], rev(y[, 3])), col = ci.col, 
                border = ci.col)
        lines(xx, y[, 1], col = line.col, lty = lty)
        if (rug) {
            Y <- x@y
            eventTimes <- Y[Y[, ncol(Y)] == 1, ncol(Y) - 1]
            rug(eventTimes, col = line.col)
        }
        return(invisible(y))
    }
    .local(x, y, ...)
}
aplot(fit,newdata=data.frame(hormon=1))

apredict <- function (object, ...) 
{
    .local <- function (object, newdata = NULL, type = c("surv", 
        "cumhaz", "hazard", "hr", "sdiff", "hdiff", "loghazard", 
        "link"), grid = FALSE, seqLength = 300, se.fit = FALSE, 
        link = NULL, exposed = incrVar(var), var, ...) 
    {
        local <- function(object, newdata = NULL, type = "surv", 
            exposed) {
            ## browser()
            tt <- object@terms
            if (is.null(newdata)) {
                X <- object@x
                XD <- object@xd
                y <- object@y
                time <- as.vector(y[, ncol(y) - 1])
            }
            else {
                lpfunc <- function(delta, fit, data, var) {
                  data[[var]] <- data[[var]] + delta
                  lpmatrix.lm(fit, data)
                }
                X <- lpmatrix.lm(object@lm, newdata)
                XD <- grad(lpfunc, 0, object@lm, newdata, object@timeVar)
                XD <- matrix(XD, nrow = nrow(X))
                if (type %in% c("hazard", "hr", "sdiff", "hdiff", 
                  "loghazard")) {
                  time <- eval(object@timeExpr, newdata)
                }
                if (object@delayed) {
                  newdata0 <- newdata
                  newdata0[[object@timeVar]] <- newdata[[object@time0Var]]
                  X0 <- lpmatrix.lm(object@lm, newdata0)
                }
                if (type %in% c("hr", "sdiff", "hdiff")) {
                  if (missing(exposed)) 
                    stop("exposed needs to be specified for type in ('hr','sdiff','hdiff')")
                  newdata2 <- exposed(newdata)
                  X2 <- lpmatrix.lm(object@lm, newdata2)
                  XD2 <- grad(lpfunc, 0, object@lm, newdata2, 
                    object@timeVar)
                  XD2 <- matrix(XD, nrow = nrow(X))
                }
            }
            beta <- coef(object)
            cumHaz = as.vector(exp(X %*% beta))
            Sigma = vcov(object)
            if (type == "link") {
                return(as.vector(X %*% beta))
            }
            if (type == "cumhaz") {
                if (object@delayed) 
                  return(cumHaz - as.vector(X0 %*% beta))
                else return(cumHaz)
            }
            if (type == "surv") {
                return(exp(-cumHaz))
            }
            if (type == "sdiff") 
                return(as.vector(exp(-exp(X2 %*% beta))) - exp(-cumHaz))
            if (type == "hazard") {
                return(as.vector(XD %*% beta) * cumHaz)
            }
            if (type == "loghazard") {
                return(as.vector(log(XD %*% beta)) + log(cumHaz))
            }
            if (type == "hdiff") {
                return(as.vector((XD2 %*% beta) * exp(X2 %*% beta) - (XD %*% 
                  beta)/time * cumHaz))
            }
            if (type == "hr") {
                cumHazRatio = exp((X2 - X) %*% beta)
                return(as.vector((XD2 %*% beta)/(XD %*% beta) * cumHazRatio))
            }
        }
        type <- match.arg(type)
        if (is.null(newdata) && type %in% c("hr", "sdiff", "hdiff")) 
            stop("Prediction using type in ('hr','sdiff','hdiff') requires newdata to be specified.")
        if (grid) {
            Y <- object@y
            event <- Y[, ncol(Y)] == 1
            time <- object@data[[object@timeVar]]
            eventTimes <- time[event]
            X <- seq(min(eventTimes), max(eventTimes), length = seqLength)[-1]
            data.x <- data.frame(X)
            names(data.x) <- object@timeVar
            newdata <- merge(newdata, data.x)
        }
        pred <- if (!se.fit) {
            local(object, newdata, type = type, exposed = exposed, 
                ...)
        }
        else {
            if (is.null(link)) 
                link <- switch(type, surv = "cloglog", cumhaz = "log", 
                  hazard = "log", hr = "log", sdiff = "I", hdiff = "I", 
                  loghazard = "I", link = "I")
            predictnl(object, local, link = link, newdata = newdata, 
                type = type, exposed = exposed, ...)
        }
        attr(pred, "newdata") <- newdata
        return(pred)
    }
    .local(object, ...)
}
environment(apredict) <- environment(stpm2)
dim(apredict(fit,newdata=data.frame(hormon=1),grid=T)) # n=300 or 299??
apredict(fit,newdata=data.frame(hormon=1),grid=T)
dim(apredict(fit,newdata=data.frame(hormon=1),grid=T,se.fit=T)) # n=300 or 299??
apredict(fit,newdata=data.frame(hormon=1),grid=T,se.fit=T)

debug(rstpm2:::numDeltaMethod)

try(suppressWarnings(detach("package:rstpm2",unload=TRUE)),silent=TRUE)
require(rstpm2)
data(brcancer)
system.time(fit2 <- stpm2(Surv(rectime/365,censrec==1)~hormon,data=brcancer,df=5))
system.time(fit3 <- pstpm2(Surv(rectime/365,censrec==1)~hormon,data=brcancer,use.gr=F))
plot(fit3,newdata=data.frame(hormon=0),type="hazard")
plot(fit2,newdata=data.frame(hormon=0),type="hazard",add=TRUE,ci=FALSE,rug=FALSE,
     line.col=2)

## penalised likelihood
brcancer$recyear <- brcancer$rectime/365
system.time(pfit1 <- pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                            logH.formula=~s(log(recyear),k=30), sp=1e-1))
system.time(fit1 <- stpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
                            logH.formula=~ns(log(recyear),df=4)))
plot(pfit1,newdata=data.frame(hormon=1))
plot(fit1,newdata=data.frame(hormon=1),lty=2,add=TRUE,ci=F)
rstpm2:::gcv(pfit1)
sps <- 10^(seq(-4,2,by=0.5))
gcvs <- sapply(sps, function(sp) {
  rstpm2:::gcv(pstpm2(Surv(recyear,censrec==1)~hormon,data=brcancer,
        logH.formula=~s(recyear,k=30), sp=sp))
})
plot(sps,gcvs,type="l",log="x")
##
system.time(fit <- rstpm2:::stpm2Old(Surv(rectime,censrec==1)~hormon,df=5,data=brcancer))
system.time(fit2 <- stpm2(Surv(rectime,censrec==1)~hormon,df=5,data=brcancer))
system.time(fit3 <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer))
##
plot(fit3,newdata=data.frame(hormon=0),type="hazard")
plot(fit2,newdata=data.frame(hormon=0),type="hazard",add=TRUE,line.col=2,ci=FALSE)
##
system.time(fit <- stpm2(Surv(rectime/365,censrec==1)~hormon,df=5,data=brcancer))
system.time(fit2 <- rstpm2:::stpm2Old(Surv(rectime/365,censrec==1)~hormon,df=5,data=brcancer))
##
system.time(fit3 <- pstpm2(Surv(rectime/365,censrec==1)~hormon,data=brcancer))
plot(fit3,newdata=data.frame(hormon=0),type="hazard")
##
plot(fit2,newdata=data.frame(hormon=0),type="hazard",add=TRUE,line.col=2,ci=FALSE)
##
plot(fit.tvc,newdata=data.frame(hormon=1),type="hr",var="hormon")
##
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                         tvc=list(hormon=3)))
anova(fit,fit.tvc) # compare with and without tvc
summary(fit.tvc <- stpm2Old(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                                      tvc=list(hormon=3)))
anova(fit,fit.tvc) # compare with and without tvc
##
plot(fit.tvc,newdata=data.frame(hormon=0),type="hr",var="hormon")
                                        # no lines method: use add=TRUE
plot(fit.tvc,newdata=data.frame(hormon=1),type="hr",var="hormon",
     add=TRUE,ci=FALSE,line.col=2)
##
## plain: identical results (good)
stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer)
stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                         logH.formula=~ns(log(rectime),3))
rstpm2:::stpm2Old(Surv(rectime,censrec==1)~hormon,data=brcancer)
## cure: identical (requires bhazard to be sensible)
rate0 <- 10^(-5+brcancer$x1/100)
(fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=2,cure=T,bhazard=rate0))
(fit2 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                         logH.formula=~nsx(log(rectime),df=2,cure=T,log=T),bhazard=rate0))
(fit3 <- rstpm2:::stpm2Old(Surv(rectime,censrec==1)~hormon,data=brcancer,cure=T,df=2,bhazard=rate0))
(fit4 <- rstpm2:::stpm2Old(Surv(rectime,censrec==1)~hormon,data=brcancer,bhazard=rate0,
                           logH.formula=~nsx(log(rectime),2,cure=T)))


##### examples #####
require(foreign)
if (FALSE) { # testing in open code
  install.packages("bbmle", repos="http://R-Forge.R-project.org")
  require(bbmle)
  brcancer=read.dta("brcancer.dta")
  brcancer=transform(brcancer,rate0=10^(-5+x1/100))
}
try(suppressWarnings(detach("package:bbmle",unload=TRUE)),silent=TRUE)

try(suppressWarnings(detach("package:rstpm2",unload=TRUE)),silent=TRUE)
## require(rstpm2)
data(brcancer)
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                     logH.formula=~nsx(log(rectime),df=3,stata=TRUE)))

brcancer <- transform(brcancer,w=10)
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                     weights=w,robust=TRUE,
                     logH.formula=~nsx(log(rectime),df=3,stata=TRUE)))


## sandwich variance estimator (from the sandwich package)

coeftest.stpm2 <- 
function (x, vcov. = NULL, df = NULL, ...) 
{
    est <- coef(x)
    if (is.null(vcov.)) 
        se <- vcov(x)
    else {
        if (is.function(vcov.)) 
            se <- vcov.(x)
        else se <- vcov.
    }
    se <- sqrt(diag(se))
    if (!is.null(names(est)) && !is.null(names(se))) {
        anames <- names(est)[names(est) %in% names(se)]
        est <- est[anames]
        se <- se[anames]
    }
    tval <- as.vector(est)/se
    pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    mthd <- "z"
    rval <- cbind(est, se, tval, pval)
    colnames(rval) <- cnames
    class(rval) <- "coeftest"
    attr(rval, "method") <- paste(mthd, "test of coefficients")
    return(rval)
}
## weights.stpm2 <- 
## function (object, ...) 
## {
##     wts <- object@weights
##     if (is.null(wts)) 
##         wts
##     else napredict(object@na.action, wts)
## }

require(sandwich)
coxph1 <- coxph(Surv(rectime,censrec==1)~hormon,data=brcancer)
update(coxph1,robust=TRUE)
sandwich(coxph1)
sandwich.stpm2(fit) # hurrah!


## require(lmtest)
## coeftest(coxph1)
## coeftest(coxph1,vcov.=sandwich(coxph1))
## coeftest(fit,sandwich(fit))


sandwich(fit)
sandwich(fit,bread.=bread.stpm2,meat.=meat.stpm2)


## some predictions
head(predict(fit,se.fit=TRUE,type="surv"))
head(predict(fit,se.fit=TRUE,type="hazard"))

## some plots
plot(fit,newdata=data.frame(hormon=0),type="hazard")

## time-varying coefficient
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                         tvc=list(hormon=3)))
anova(fit,fit.tvc) # compare with and without tvc

plot(fit.tvc,newdata=data.frame(hormon=0),type="hr",var="hormon")
                                        # no lines method: use add=TRUE
plot(fit.tvc,newdata=data.frame(hormon=1),type="hr",var="hormon",
     add=TRUE,ci=FALSE,line.col=2)

plot(fit.tvc,newdata=data.frame(hormon=0),type="sdiff",var="hormon")

plot(fit.tvc,newdata=data.frame(hormon=0),type="hdiff",var="hormon")

plot(fit.tvc,newdata=data.frame(hormon=0),type="hazard")
plot(fit.tvc,newdata=data.frame(hormon=1),type="hazard",line.col=2,ci=FALSE,add=TRUE)
## trace("predict", browser, exit=browser, signature = "stpm2")

set.seed(10101)
brcancer <- transform(brcancer, x=rlnorm(nrow(brcancer)))
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                     tvc.formula=~hormon:nsx(log(rectime),df=3)))



## cure model
## cf. http://www.pauldickman.com/survival/solutions/q37.do
### Data setup
require(foreign)
colon <- read.dta("http://www.pauldickman.com/survival/colon.dta")
popmort <- read.dta("http://www.pauldickman.com/survival/popmort.dta")
brcancer <- read.dta("http://www.stata-press.com/data/r11/brcancer.dta")
popmort <- transform(popmort, age=`_age`, year=`_year`, `_age`=NULL, `_year`=NULL)

save(colon,file="c:/usr/src/R/rstpm2/pkg/data/colon.rda")
save(popmort,file="c:/usr/src/R/rstpm2/pkg/data/popmort.rda")
save(brcancer,file="c:/usr/src/R/rstpm2/pkg/data/brcancer.rda")

## require(rstpm2)
popmort2 <- transform(popmort,exitage=age,exityear=year,age=NULL,year=NULL)
colon2 <- within(colon, {
  status <- ifelse(surv_mm>120.5,1,status)
  tm <- pmin(surv_mm,120.5)/12
  exit <- dx+tm*365.25
  sex <- as.numeric(sex)
  exitage <- pmin(floor(age+tm),99)
  exityear <- floor(yydx+tm)
})
colon2 <- merge(colon2,popmort2)

## compare relative survival without and with cure 
summary(fit0 <- stpm2(Surv(tm,status %in% 2:3)~I(year8594=="Diagnosed 85-94"),
                     data=colon2,
                     bhazard=colon2$rate, df=5)) ## CHECKED: same year8594 estimate as Stata
head(predict(fit0))
## estimate of failure at the end of follow-up
1-predict(fit0,data.frame(year8594 = unique(colon2$year8594),tm=max(colon2$tm)),type="surv",se.fit=TRUE)
plot(fit0,newdata=data.frame(year8594 = "Diagnosed 85-94"),ylim=0:1)
plot(fit0,newdata=data.frame(year8594 = "Diagnosed 75-84"),add=TRUE,line.col="red",rug=FALSE)
##
summary(fit <- stpm2(Surv(tm,status %in% 2:3)~I(year8594=="Diagnosed 85-94"),
                     data=colon2,
                     bhazard=colon2$rate,
                     df=5,cure=TRUE))
head(predict(fit))
## cure fractions (I need to add this to the predict function)
1-predict(fit,data.frame(year8594 = unique(colon2$year8594),tm=max(colon2$tm)),type="surv",se.fit=TRUE)
newdata1 <- data.frame(year8594 = "Diagnosed 85-94")
plot(fit,newdata=newdata1,add=TRUE,ci=FALSE,lty=2,rug=FALSE)
plot(fit,newdata=data.frame(year8594="Diagnosed 75-84"),add=TRUE,rug=FALSE,line.col="red",ci=FALSE,lty=2)

plot(fit,newdata=newdata1,type="hazard")
plot(fit,newdata=newdata1,type="cumhaz")


## http://www.pauldickman.com/survival/r/melanoma.relsurv.r
library(foreign)
library(survival)
library(relsurv)
# Download rates files from http://www.mortality.org/
# # 6. Life Tables By year of death (period) 1x1
# Save tables by gender in text files
# The transrate.hmd command translate these to R ratetables
Finlandpop <- transrate.hmd("c:/usr/tmp/mltper_1x1.txt","c:/usr/tmp/fltper_1x1.txt")

## The relsurv package requires time in days (exit and dx are dates of exit and diagnosis)
colon3 <- transform(colon2,tm.dd=as.numeric(exit-dx))
colon3$sex <- ifelse(colon2$sex==1,"male","female")
as.date <- function(x)
  if (inherits(x,"Date")) as.date(as.numeric(x)+3653) else date::as.date(x)
model1 <- rs.surv(Surv(tm.dd,status %in% 2:3)~year8594+ratetable(age=(X_age+0.5)*365.25,sex=sex,year=as.date(exit)),colon3,ratetable=Finlandpop)
plot(model1,lty=1:2)



                 
oldx <- 0:100
oldy <- (oldx-50)^2
oldy[c(20,30)] <- 0
old <- data.frame(x=oldx,y=oldy)
predict(lm(y~nsx(x,knots=c(25,50,75,95)),old)) # as per Stata
newx <- seq(min(oldx)/1.05,max(oldx)*1.05,length=101)
new <- data.frame(x=newx)
plot(oldx,oldy)
predict(lm(y~nsx(x,df=5,cure=TRUE),old))
sum(oldy)
terms(lm(y~nsx(x,df=5,cure=TRUE),old))
lm(y~nsx(x,df=5),old)


lines(newx,
      predict(lm(y~nsx(x,df=4,cure=FALSE),old),newdata=new),
      type="l") # oops
lines(newx,
      predict(lm(y~nsx(x,df=3),old),newdata=new),
      lty=2)


summary(fit <- stpm2(Surv(tm,status %in% 2:3)~I(year8594=="Diagnosed 85-94"),
                     data=colon2,
                     bhazard=colon2$rate,
                     logH.formula=~nsx(log(tm),df=6,stata=TRUE))) # okay
summary(fit <- stpm2(Surv(tm,status %in% 2:3)~I(year8594=="Diagnosed 85-94"),
                     data=colon2,
                     logH.formula=~nsx(log(tm),df=6,stata=TRUE))) # okay

## Stata
## stata.knots=c(4.276666164398193, 6.214608192443848, 6.7833251953125, 7.806289196014404)
stataKnots <- function(x,df) {
  intKnots <- round((1:(df-1))/df,2) # yes, Paul implicitly rounded to 2 dp
  logx <- log(x)
  c(min(logx),quantile(logx,intKnots,type=2),max(logx))
}
stata.knots <- stataKnots(subset(brcancer,censrec==1)$rectime,3)
## sapply(1:9,function(type) log(quantile(subset(brcancer,censrec==1)$rectime,c(0.33,0.67),type=type)))
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                     logH.args=list(knots=stata.knots[2:3],
                       Boundary.knots=stata.knots[c(1,4)])))
## formula specification for logH
summary(stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
              logH.formula=~ns(log(rectime),df=3)))

pred <- predict(fit.tvc,newdata=data.frame(hormon=0:3),grid=TRUE,se.fit=TRUE,type="cumhaz")
pred.all <- cbind(pred,attr(pred,"newdata"))
require(lattice)
xyplot(Estimate ~ rectime, data=pred.all, group=hormon,type="l",xlab="Time")


## relative survival
brcancer <- transform(brcancer,rate0=10^(-5+x1/100))
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,bhazard=brcancer$rate0,df=3))
head(predict(fit,se.fit=TRUE))

## delayed entry
brcancer2 <- transform(brcancer,startTime=ifelse(hormon==0,rectime*0.5,0))
## debug(stpm2)
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
                     logH.formula=~nsx(log(rectime),df=3,stata=TRUE)))
head(predict(fit,se.fit=TRUE))
## delayed entry and tvc
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
                     tvc.formula=~hormon:nsx(log(rectime),df=3,stata=TRUE)))
head(predict(fit,se.fit=TRUE))



## multiple time scales
brcancer <- transform(brcancer,recyr=rectime/365.25)
## predictions from a simple model
summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon+x1,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50))))
head(predict(fit))
grid.x1 <- with(brcancer, seq(40,70,length=300))
newdata0 <- with(brcancer, data.frame(recyr=5,x1=grid.x1,hormon=0))
matplot(grid.x1,
        predict(fit,type="hr",newdata=newdata0,var="hormon",se.fit=TRUE), type="l")
## predictions with multiple time scales
summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50)),
                     tvc.formula=~hormon:nsx(log(recyr+x1),df=2)))
matplot(grid.x1,
        predict(fit,type="hr",newdata=newdata0,var="hormon",se.fit=TRUE), type="l")





brcancer <- transform(brcancer,recyr=rectime/365.25,entry=recyr/2)
summary(fit <- stpm2(Surv(entry,recyr,censrec==1)~hormon,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50)),
                     tvc.formula=~hormon:nsx(log(recyr+x1),df=2)))


summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon+x1,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50))))


plot(grid.x1,
     predict(fit,type="hr",newdata=newdata0,var="hormon",se.fit=TRUE)$fit, type="l")

plot(fit,newdata=data.frame(hormon=0,x1=50),var="hormon",type="hr")

head(predict(fit,type="hazard",newdata=newdata0))
head(predict(fit,type="hazard",newdata=transform(newdata0,hormon=1)))



newdata0 <- with(brcancer, data.frame(recyr=5+1,x1=grid.x1-1,hormon=0))
predict(fit,type="hr",newdata=newdata0,var="hormon")

summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon+x1,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50)),tvc=list(hormon=3)))

brcancer <- transform(brcancer, startAge=x1, endAge=x1+rectime/365)
summary(fit <- stpm2(Surv(startAge,endAge,censrec==1)~hormon,data=brcancer,
                     logH.formula=~nsx(log(endAge),df=3,centre=log(50)),tvc=list(hormon=3)))


## some simulated data: H_{weibull}(t)=(t/b)^a
n <- 1000
sim1 <- data.frame(age=seq(20,70,length=n),x=rep(0:1,each=n/2))
y <- rweibull(1000,shape=1,scale=1)



with(brcancer, plot(density(x1[censrec==1])))

summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon,data=brcancer,logH.formula=~nsx(log(recyr),df=3,stata=TRUE)))




brcancer <- transform(brcancer,ageStart=rnorm(length(rectime),50,5))
brcancer <- transform(brcancer,ageStop=ageStart+rectime)
summary(fit <- stpm2(Surv(ageStart,ageStop,censrec==1)~hormon,data=brcancer,df=3))

brcancer3 <- transform(brcancer,startTime=ifelse(censrec==1,0,10))
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=subset(brcancer,rectime>10),df=3))
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=subset(brcancer3,rectime>10),df=3))

## check the performance time
refresh
require(rstpm2)
data(brcancer)
brcancer10 = do.call("rbind",lapply(1:10,function(i) brcancer))
system.time(summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,df=3,data=brcancer10)))
system.time(summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,df=3,data=brcancer10)))
system.time(summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer10, logH.formula=~ns(log(rectime),df=4)+hormon:ns(log(rectime),df=3))))
system.time(summary(fit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer10)))
system.time(summary(fit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer10, smooth.formula=~s(log(rectime))+s(log(rectime),by=hormon))))

fit <- pstpm2(Surv(rectime,censrec==1)~hormon,data=brcancer10,trace=1)

fit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer10, smooth.formula=~s(log(rectime))+s(log(rectime),by=hormon),trace=1,sp.init=c(1,1), reltol=list(outer=1e-5,search=1e-10,final=1e-10))

system.time(summary(fit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer10, 
                                  smooth.formula=~s(log(rectime))+s(log(rectime),by=hormon),
                                  sp=c(0.006,0.0031),trace=1,outer_optim=2,criterion="GCV",
                                  reltol=list(outer=1e-5,search=1e-10,final=1e-10))))
## > fit@sp
## [1] 0.06104312 0.31430954

system.time(fit <- pstpm2(Surv(rectime,censrec==1)~1,data=brcancer10, smooth.formula=~s(log(rectime))+s(log(rectime),by=hormon),sp=c(1,1)))


nsx(1:10,df=3) - ns(1:10,df=3)
nsx(1:10,df=3,centre=3)
nsx(1:10,df=3,centre=3,Boundary.knots=c(2,8),derivs=c(1,1))
nsx(1:10,df=3,cure=TRUE)
nsxDeriv(1:10,df=3) - nsDeriv(1:10,df=3)
nsxDeriv(1:10,df=3,centre=5,derivs=c(1,1))
nsxDeriv(1:10,df=3,centre=5,cure=TRUE)

nsDeriv(1:10,df=3) - nsDeriv2(1:10,df=3)

## bug with calling mle2
require(bbmle)
mle2a <- function(...)
  mle2(...)
## some data
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)
## some fits
(fit0 <- mle2(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d)) # okay
(fit0.2 <- mle2(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d,
              control=list(parscale=2))) # okay
(fit1 <- mle2a(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d)) # okay
(fit1.2 <- mle2a(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d,
              control=list(parscale=2))) # FAILS

## end of examples ##




## ## * stata
## cd c:\Users\marcle\Documents\Home\
## clear
## webuse brcancer
## use brcancer
## stset rectime, f(censrec==1)
## cap program drop dopredictions
## program define dopredictions
##   preserve
##   predict hr, hrnumerator(hormon 1) ci
##   predict haz, hazard ci
##   predict surv, surv ci
##   predict sdiff, sdiff1(hormon 1) ci
##   list hr* in 1/5
##   list haz* surv* in 1/5
##   list sdiff* in 1/5
##   restore
## end

## * basic model
## stpm2 hormon, df(3) scale(h)
## dopredictions

## * cure 
## gen rate0=10^(-5+x1/100)
## stpm2 hormon, df(3) scale(h) cure bhazard(rate0)
## dopredictions

## * tvc
## stpm2 hormon, df(3) tvc(hormon) dftvc(3) scale(h)
## dopredictions

## * delayed entry
## preserve
##   replace _t0 = rectime*0.5 if hormon==0
##   stpm2 hormon, df(3) scale(h)
##   dopredictions
## restore

## * relative survival
## preserve  
##   gen rate0=10^(-5+x1/100)
##   stpm2 hormon, df(3) scale(h) bhazard(rate0)
##   dopredictions
## restore

## * test speed
## clear all
## set mem 100m
## use brcancer
## stset rectime, f(censrec==1)
## expand 100
## timer clear
## timer on 1
## stpm2 hormon, df(3) scale(h)
## timer off 1
## timer list
 

## hazard.pm = function(obj,tm,X,XD) # obj$par
## {
##   Xlocal=predict(X,newx=log(tm))
##   XDlocal=predict(XD,newx=log(tm))
##   with(obj,
##        c((XDlocal %*% par)/tm*exp(Xlocal %*% par)))
## }
## with(list(df=df,x=seq(0,3,length=100)[-1]),
##      {
##        plot(x,hazard.pm(fit,x,X,XD),type="l",ylim=c(0,2))
##        lines(x,dweibull(x,shape=1)/pweibull(x,shape=1,lower=FALSE),lty=2)
##      })
## ##
## require(deSolve)
## temp <- as.data.frame(ode(y=0,times=seq(0,10,length=100)[-1],
##                           func=function(t,state,parameters=NULL) list(exp(sin(2*pi*log(t))))))
## plot(temp,type="l")
## temp <- transform(temp, cum=`1`,logcum=log(`1`))
## with(temp,plot(log(time),logcum))
## temp1 <- temp[-1,]
## fit <- glm(log(cum)~log(time)+sin(2*pi*log(time))+cos(2*pi*log(time)),data=temp1)
## lines(log(temp1$time),predict(fit))
## ## In summary:
## ## we can model using sine and cosine terms for the log-cumulative hazard - for log(time).
