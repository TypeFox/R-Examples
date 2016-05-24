## This example is provided to support the comparison in the
## discussion section of the supporting Gramacy & Polson paper

## load the data
data(pima)
X <- as.matrix(pima[,-9])
y <- as.numeric(pima[,9])
     
## pre-normalize to match the comparison in the paper
one <- rep(1, nrow(X))
normx <- sqrt(drop(one %*% (X^2)))
X <- scale(X, FALSE, normx)

## compare to the GLM fit
fit.logit <- glm(y~X, family=binomial(link="logit"))

## do the Gibbs sampling
T <- 10000
system.time(out6 <- reglogit(T, y, X, nu=0, nup=NULL, normalize=FALSE))
     
## plot the posterior distribution of the coefficients
burnin <- (1:(T/10)) 
boxplot(out6$beta[-burnin,], ylab="posterior",
        xlab="coefficients", bty="n",  names=c("mu", paste("b", 1:8, sep="")))
abline(h=0, lty=2)

## add in GLM fit and MAP with legend
points(fit.logit$coef, col=2, pch=17)
points(out6$map$beta, pch=19, col=3)
legend("topright", c("MLE", "MAP"), col=2:3, pch=c(17,19))

## using dRUM
system.time(fit.dRUM <- gibbs.dRUM(T, cbind(rep(1,length(y)), X), y))

## add to boxplot
boxplot(fit.dRUM$beta[-burnin,], add=TRUE, col=3, names=rep("",9))

## calculate effective sizes
ess.dRUM <- apply(fit.dRUM$beta[-burnin,], 2, effectiveSize)
ess <- apply(out6$beta[-burnin,], 2, effectiveSize)
mean(ess.dRUM/ess)
## ~ 20%

## These are the times (but they need to be re-run when Arnor is less busy
## > 144/753
## [1] 0.1912351

