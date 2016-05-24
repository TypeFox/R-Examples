library(coxme)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

#
# Check out the code on a simulated data set, 
#  arose out a concern that the first derivatives were wrong
#  when the model had random slopes and intercepts
# (The code was ok, but another test case is always good.)
# The heart of the test is to re-compute  quantities using
#   a combination of coxph or matrix algebra, and check these with 
#   the coxme answers.
#
approx <- "efron"
set.seed(9291978)
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
fixf <- coxph(Surv(futime, status) ~ factor(site)*trt1, trdata,
              method=approx)

fit1 <- coxme(Surv(futime, status) ~ trt1 + (1 + trt1 | site), trdata,
              ties=approx)
fit2 <- coxme(Surv(futime, status) ~ trt2 + (1 + trt2 | site), trdata,
              ties=approx)
# Fit1 and fit2 are actually different models, with different 
#  iteration paths and maxima.

# Test out loglik
eta <- fit1$linear
fit1b <- coxph(Surv(futime, status) ~ eta, trdata, init=1, iter=0,
               method=approx)
aeq(fit1b$loglik[2], fit1$loglik[3]+fit1$penal)

# First derivatives
temp <- VarCorr(fit1)[[1]]
sigma <- diag(c(rep(temp[1],4), rep(temp[4],4)))
sigma[cbind(1:4,5:8)] <- temp[3]* sqrt(temp[1] * temp[4])
sigma[cbind(5:8,1:4)] <- temp[3]* sqrt(temp[1] * temp[4])
pen <- matrix(0., 9,9)
pen[1:8, 1:8] <- solve(sigma)
pcoef <- c(unlist(ranef(fit1)), fixef(fit1))
aeq(pcoef %*% pen %*% pcoef/2, fit1$penal)

xx <- 1* with(trdata, cbind(site==1, site==2, site==3, site==4))
xx <- cbind(xx, xx*trdata$trt1)

cfit1 <- coxph(Surv(futime, status) ~ xx + trt1, iter=0,
               method=approx, init=pcoef, trdata)
aeq(cfit1$log, fit1b$log) #check my typing

# Using coxph.detail, I find the first and second derivatives
#  via a completely different code path than coxme
dt1 <- coxph.detail(cfit1)
first <- colSums(dt1$score) - pcoef%*%pen
aeq(first, fit1$u)  #do they agree?

i2 <- apply(dt1$imat, 1:2, sum) + pen  #information is -1*second deriv
aeq(as.matrix(fit1$hmat, one=F), as.matrix(gchol(i2), ones=F))
