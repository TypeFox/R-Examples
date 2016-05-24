## This example is provided to support the comparison in the
## discussion section of the supporting Gramacy & Polson paper

## get the data
spam <- read.table("spamdata.txt")

## extract X and y -- using logarithm'd predictors
X <- as.matrix(spam[,-58])
y <- drop(spam[,58])

## compare to the GLM fit
fit.logit <- glm(y~X, family=binomial(link="logit"))

## do the Gibbs sampling
T <- 10000
system.time(out6 <- reglogit(T, y, X, nu=0, nup=NULL,normalize=FALSE))
     
## using dRUM
system.time(fit.dRUM <- gibbs.dRUM(T, cbind(rep(1,length(y)), X), y))

burnin <- (1:(8*T/10)) 

## calculate effective sizes
ess.dRUM <- apply(fit.dRUM$beta[-burnin,], 2, effectiveSize)
ess <- apply(out6$beta[-burnin,], 2, effectiveSize)
mean(ess.dRUM/ess)
## 7.6%

## time comparison
## > 343/1833
## [1] 0.1871249

## rejection rate
length(unique(fit.dRUM$beta[,1]))/T
## 0.007
