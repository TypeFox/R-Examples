###################################################
### chunk number 1: 
###################################################
#line 6 "cap5.Rnw"
options(prompt="R> ")
options(width=80)


###################################################
### chunk number 2: simGBM
###################################################
#line 138 "cap5.Rnw"
require(sde)
set.seed(123)
sigma <- 0.5
mu <- 1
S <- sde.sim(X0=100, model="BS", theta=c(mu, sigma),T=100,N=10000)


###################################################
### chunk number 3: 
###################################################
#line 146 "cap5.Rnw"
X <- diff(log(S))
sigma.hat <- sqrt( var(X)/deltat(S) )
alpha.hat <- mean(X)/deltat(S)
mu.hat <- alpha.hat + 0.5 * sigma.hat^2
sigma.hat
mu.hat


###################################################
### chunk number 4: qmle eval=FALSE
###################################################
## #line 245 "cap5.Rnw"
## require(yuima)
## diff.matrix <- matrix(c("alpha1","alpha2", "1", "1"), 2, 2)
## 
## drift.c <- c("-1*beta1*x1", "-1*beta2*x2", "-1*beta2", "-1*beta1")
## drift.matrix <- matrix(drift.c, 2, 2)
## 
## ymodel <- setModel(drift=drift.matrix, diffusion=diff.matrix, time.variable="t",
##                    state.variable=c("x1", "x2"), solve.variable=c("x1", "x2"))
## n <- 100
## ysamp <- setSampling(Terminal=(n)^(1/3), n=n)
## yuima <- setYuima(model=ymodel, sampling=ysamp)
## 
## set.seed(123)
## truep <- list(alpha1=0.5, alpha2=0.3, beta1=0.6, beta2=0.2)
## yuima <- simulate(yuima, xinit=c(1, 1), true.parameter=truep)
## 
## opt <- qmle(yuima, start=list(alpha1=0.7, alpha2=0.2, beta1=0.8, beta2=.2))


###################################################
### chunk number 5: 
###################################################
#line 265 "cap5.Rnw"
#line 245 "cap5.Rnw#from line#265#"
require(yuima)
diff.matrix <- matrix(c("alpha1","alpha2", "1", "1"), 2, 2)

drift.c <- c("-1*beta1*x1", "-1*beta2*x2", "-1*beta2", "-1*beta1")
drift.matrix <- matrix(drift.c, 2, 2)

ymodel <- setModel(drift=drift.matrix, diffusion=diff.matrix, time.variable="t",
                   state.variable=c("x1", "x2"), solve.variable=c("x1", "x2"))
n <- 100
ysamp <- setSampling(Terminal=(n)^(1/3), n=n)
yuima <- setYuima(model=ymodel, sampling=ysamp)

set.seed(123)
truep <- list(alpha1=0.5, alpha2=0.3, beta1=0.6, beta2=0.2)
yuima <- simulate(yuima, xinit=c(1, 1), true.parameter=truep)

opt <- qmle(yuima, start=list(alpha1=0.7, alpha2=0.2, beta1=0.8, beta2=.2))
#line 266 "cap5.Rnw"


###################################################
### chunk number 6: 
###################################################
#line 268 "cap5.Rnw"
opt@coef
unlist(truep)


###################################################
### chunk number 7: interest eval=FALSE
###################################################
## #line 277 "cap5.Rnw"
## library(Ecdat)
## library(sde)
## data(Irates)
## X <- Irates[,"r1"]
## plot(X)


###################################################
### chunk number 8: 
###################################################
#line 286 "cap5.Rnw"
par(mar=c(3,3,1,1))
#line 277 "cap5.Rnw#from line#287#"
library(Ecdat)
library(sde)
data(Irates)
X <- Irates[,"r1"]
plot(X)
#line 288 "cap5.Rnw"


###################################################
### chunk number 9: 
###################################################
#line 355 "cap5.Rnw"
sde.sim(model="CIR", theta=c(1, 0.3, .1))


###################################################
### chunk number 10: 
###################################################
#line 359 "cap5.Rnw"
CIR.loglik <- function(theta1,theta2,theta3) {
 n <- length(X)
 dt <- deltat(X)
 -sum(dcCIR(x=X[-1], Dt=dt, x0=X[-n], theta=c(theta1,theta2,theta3),
   log=TRUE))
}


###################################################
### chunk number 11: mleCIR
###################################################
#line 368 "cap5.Rnw"
fit <-  mle(CIR.loglik, start=list(theta1=.1,  theta2=.1,theta3=.3),  method="L-BFGS-B",lower=rep(1e-3,3), upper=rep(1,3)) 
coef(fit)


###################################################
### chunk number 12: 
###################################################
#line 413 "cap5.Rnw"
library(Ecdat)
data(Irates)
X <- Irates[,"r1"]

Y <- X[-1]
stage1 <- lm(diff(X) ~ I(1/Y) + Y + I(Y^2))
coef(stage1)

eps2 <- residuals(stage1)^2
mod <- nls(eps2 ~ b0 + b1*Y +b2*Y^b3, start=list(b0=1, b1=1, b2=1, b3=.5),lower=rep(1e-5,4), upper=rep(2,4), algorithm="port")

w <- predict(mod)
stage2 <- lm(diff(X) ~ I(1/Y) + Y + I(Y^2), weights=1/w)

coef(stage2)
coef(mod)


###################################################
### chunk number 13: 
###################################################
#line 432 "cap5.Rnw"
summary(stage2)
summary(mod)


###################################################
### chunk number 14: 
###################################################
#line 450 "cap5.Rnw"
require(fImport)
data <- yahooSeries("ENI.MI", from="2009-01-01", to="2009-12-31")
S <- data[, "ENI.MI.Close"]
X <- returns(S)


###################################################
### chunk number 15: ENI
###################################################
#line 458 "cap5.Rnw"
require(quantmod)
lineChart(S,layout=NULL,theme="white") 
lineChart(X,layout=NULL,theme="white")


###################################################
### chunk number 16: 
###################################################
#line 465 "cap5.Rnw"
par(mfrow=c(2,1))
#line 458 "cap5.Rnw#from line#466#"
require(quantmod)
lineChart(S,layout=NULL,theme="white") 
lineChart(X,layout=NULL,theme="white")
#line 467 "cap5.Rnw"


###################################################
### chunk number 17: levyfit eval=FALSE
###################################################
## #line 474 "cap5.Rnw"
## library(fBasics)
## nFit(X)
## nigFit(X,trace=FALSE)
## hypFit(X,trace=FALSE)
## ghFit(X,trace=FALSE)


###################################################
### chunk number 18: 
###################################################
#line 483 "cap5.Rnw"
par(mfrow=c(2,2))
par(mar=c(3,3,3,1))
grid <- NULL
#line 474 "cap5.Rnw#from line#486#"
library(fBasics)
nFit(X)
nigFit(X,trace=FALSE)
hypFit(X,trace=FALSE)
ghFit(X,trace=FALSE)
#line 487 "cap5.Rnw"


###################################################
### chunk number 19: 
###################################################
#line 948 "cap5.Rnw"
library(sde)
AhnGao.loglik <- function(theta1,theta2,theta3) {
 n <- length(X)
 dt <- deltat(X)
 -sum(dcCIR(x=1/X[-1], Dt=dt, x0=1/X[-n], theta=c(theta1,theta2,theta3),
   log=TRUE))+2*sum(log(X[-1]))
}


###################################################
### chunk number 20: 
###################################################
#line 958 "cap5.Rnw"
library(Ecdat)
data(Irates)
X <- Irates[,"r1"]
fit <-  mle(AhnGao.loglik, start=list(theta1=.1,  theta2=.1,theta3=.3),  method="L-BFGS-B",lower=rep(1e-3,3), upper=rep(1,3)) 


###################################################
### chunk number 21: 
###################################################
#line 965 "cap5.Rnw"
coef(fit)


