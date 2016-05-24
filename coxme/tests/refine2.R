#
# Check out the code on a simulated data set.  On first cut it appeared to
#  do very badly, this turned out to be a bug in coxfit6d.c, when there
#  were random slopes it miscalculated the linear predictor.
#
library(coxme)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

set.seed(1979)
mkdata <- function(n, beta=c(.4, .1), sitehaz=c(.5,1.5, 2,1)) {
    nsite <- length(sitehaz)
    site <- rep(1:nsite, each=n)
    trt1 <- rep(0:1, length=n*nsite)
    hazard <- sitehaz[site] + beta[1]*trt1 + beta[2]*trt1 * (site-mean(site))
    stime <- rexp(n*nsite, exp(hazard))
    q80 <- quantile(stime, .8)
    data.frame(site=site,
               trt1 = trt1,
               trt2 = 1-trt1,
               futime= pmin(stime, q80),
               status= ifelse(stime>q80, 0, 1),
               hazard=hazard
               )
    }

trdata <- mkdata(150)  #150 enrolled per site
#fixf <- coxph(Surv(futime, status) ~ factor(site)*trt1, trdata)
#pfit <- pyears(Surv(futime, status) ~ site + trt1, trdata)
#(pfit$event/sum(pfit$event))/ (pfit$pyears/sum(pfit$pyears))

set.seed(50)
#nsim <- 500
nsim <- 20  # speedup for CRAN

fit <- coxme(Surv(futime, status) ~ trt2 + (1 + trt2 | site), trdata,
             refine.n=nsim, refine.detail=TRUE)
debug <- fit$refine.detail

# Recreate the variance-covariance and sim data that was used
set.seed(50)
bmat <- matrix(rnorm(8*nsim), nrow=8)
hmatb <- fit$hmat[1:8, 1:8]
hmat.inv <- as.matrix(solve(hmatb))

chidf <- coxme.control()$refine.df

b2 <- backsolve(hmatb, bmat, upper=TRUE)
htran <- as(hmatb, "dtCMatrix")  #verify that backsolve works correctly
all.equal(as.matrix(htran %*% b2), bmat, check.attributes=FALSE)

b2 <- scale(b2, center=F, scale=sqrt(rchisq(nsim, df=chidf)/chidf))
b3 <- b2 + unlist(ranef(fit)) 
aeq(b3, debug$bmat)


temp <- VarCorr(fit)[[1]]
sigma <- diag(c(rep(temp[1],4), rep(temp[4],4)))
sigma[cbind(1:4,5:8)] <- temp[3]* sqrt(temp[1] * temp[4])
sigma[cbind(5:8,1:4)] <- temp[3]* sqrt(temp[1] * temp[4])

if (!is.null(debug))
    all.equal(as.matrix(gchol(sigma), ones=F), as.matrix(debug$gkmat, ones=F))

coxll <- double(nsim)
for (i in 1:nsim) {
    lp <- trdata$trt2 * fixef(fit) + b3[trdata$site,i] +
           b3[trdata$site +4,i] * trdata$trt2
    tfit <- coxph(Surv(futime, status) ~ offset(lp), trdata)
    coxll[i] <- tfit$loglik
    }
if (!is.null(debug)) all.equal(coxll, debug$log)

# How does the direct average compare to the IPL?
# (This is not guarranteed to be accurate, just curious)
constant <- .5*(log(2*pi) + sum(log(diag(gchol(sigma)))))
fit$log[2] + c(IPL=0, sim=log(mean(exp(coxll-fit$log[2]))) - constant)

# Compute the Taylor series for the IPL
bhat <- unlist(ranef(fit))
b.sig <- t(b2) %*% fit$hmat[1:8, 1:8]  # b times sqrt(H)
taylor <- rowSums(b.sig^2)/2   # vector of b'H b/2
aeq(taylor, debug$penalty2)

# Look at how well the Taylor series does
pen <- solve(sigma)
simpen <- colSums(b3 *(pen%*%b3))/2 #penalty for each simulated b
aeq(simpen, debug$penalty1)

#plot(coxll- simpen, fit$log[3] - taylor, 
#     xlab="Actual PPL for simulated points", 
#     ylab="Taylor series approximation")
#abline(0,1)
#

#And now I need the Gaussian integration contant and the t-density
require(mvtnorm)
tdens <- dmvt(t(b3), delta=unlist(ranef(fit)), sigma=hmat.inv, df=chidf)
aeq(tdens, debug$tdens)
 
# The normalization constant for the Gaussian penalty
#  
temp <- determinant(sigma, log=TRUE)
gnorm <- -(4*log(2*pi) + .5*temp$sign *temp$modulus)
aeq(gnorm, debug$gdens)

m2 <- fit$loglik[2]
errhat <- exp(coxll+gnorm -(simpen + tdens +m2)) - 
         exp(fit$log[3]+ gnorm -(taylor + tdens + m2))
if (!is.null(debug)) {
    aeq(errhat, debug$errhat)
    }
