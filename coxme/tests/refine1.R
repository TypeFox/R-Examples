library(coxme)
# Test of refine.n option
#  Matches the first example in the laplace vignette
#
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)
nsim <- 100
set.seed(1953)  #an auspicious birth year :-)
mkdata <- function(n, beta=c(.4, .1), sitehaz=c(.5,1.5, 2,1)) {
    nsite <- length(sitehaz)
    site <- rep(1:nsite, each=n)
    trt1 <- rep(0:1, length=n*nsite)
    hazard <- sitehaz[site] + beta[1]*trt1 + beta[2]*trt1 * (site-mean(site))
    stime <- rexp(n*nsite, exp(hazard))
    q80 <- quantile(stime, .8)
    data.frame(site=site,
               trt = trt1,
               trt2 = 1-trt1,
               futime= pmin(stime, q80),
               status= ifelse(stime>q80, 0, 1),
               hazard=hazard
               )
}
trdata <- mkdata(150)  #150 enrolled per site

set.seed(12345)

fit1 <- coxme(Surv(futime, status) ~  trt + (1| site/trt), trdata,
              refine.n=nsim, refine.detail=TRUE)

hmatb <- fit1$hmat[1:12, 1:12]
hmat.inv <- as.matrix(solve(hmatb))

chidf <- coxme.control()$refine.df
set.seed(12345)
bmat <- matrix(rnorm(12*nsim), ncol=nsim)

b2 <- backsolve(hmatb, bmat, upper=TRUE)
htran <- as(hmatb, "dtCMatrix")  #verify that backsolve works correctly
all.equal(as.matrix(htran %*% b2), bmat, check.attributes=FALSE)

b2 <- scale(b2, center=F, scale=sqrt(rchisq(nsim, df=chidf)/chidf))
b3 <- b2 + unlist(ranef(fit1)) 
aeq(b3, fit1$refine.detail$bmat)

# Check the fitted Cox models
logvec <- double(nsim)
indx <- -1 + trdata$trt + 2*trdata$site  #random treatment effect index
for (i in 1:nsim) {
    eta <- fixef(fit1) * trdata$trt + b3[indx, i] +
            b3[trdata$site+8,i] 
    tfit <- coxph(Surv(futime, status) ~ offset(eta), data= trdata)
    logvec[i] <- tfit$loglik[1] 
}
aeq(logvec, fit1$refine.detail$loglik)

# The penalties
temp <- unlist(VarCorr(fit1))
Ainv <- diag(rep(1/temp, c(8,4)))
p1 <- diag(t(b3) %*% Ainv %*% b3)/2
aeq(p1, fit1$refine.detail$penalty1)

p2 <- mahalanobis(t(b3), center=unlist(ranef(fit1)), cov=hmat.inv)
p3 <- mahalanobis(t(b2), center=rep(0,12), cov=hmat.inv)
aeq(p2, p3)
aeq(p2/2, fit1$refine.detail$penalty2)

# The log-density of the t
require(mvtnorm)
tdens <- dmvt(t(b3), delta=unlist(ranef(fit1)), sigma=hmat.inv, df=chidf)
aeq(tdens, fit1$refine.detail$tdens)
 
# The normalization constant for the Gaussian penalty
#  Determinants are easy for a diagonal matrix
gnorm <- -(6*log(2*pi) - .5*sum(log(diag(Ainv))))

# My errhat vector
t1 <- logvec + gnorm - (tdens + p1 + fit1$loglik[2])
t2 <- fit1$loglik[3] + gnorm - (tdens + p2/2 + fit1$loglik[2])
errhat <- exp(t1) - exp(t2)
aeq(errhat, fit1$refine.detail$errhat)
