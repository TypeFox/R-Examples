library(coxme)
options(na.action='na.exclude', contrasts=c('contr.treatment', 'contr.poly'))
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

#
# Same data set as slope1
#
set.seed(56)
n.subject <- seq(180, by=21, length=9) # number of subjects
slope <- sort(-.5 + rnorm(9, sd=.5))         # true treament effects

inst <- rep(1:9, n.subject)
n <- length(inst)
simdata <- data.frame(id=1:n, inst=inst,
                      trt= rep(0:1, length=n),
                      age= runif(n, 40, 70))
#risk goes up 30%/decade of age
simdata$hazard <- .8* exp(simdata$trt * rep(slope, n.subject) +
                          (simdata$age-55) * .03)

rtime <- function(hazard, censor=c(1,2)) {
    stime <- rexp(length(hazard), rate=hazard)
    ctime <- runif(length(hazard), censor[1], censor[2])
    list(time= pmin(stime, ctime), status=1*(stime <=ctime))
    }
temp <- rtime(simdata$hazard)
simdata$time <- temp$time
simdata$status <- temp$status

#
# Test out the refine.n code, using the simdata
#  A simple diagonal variance
#
# For original testing we had nsim=100, changed to 10 for a CRAN speedup
nsim <- 10
var  <- .3   #sizeable

set.seed(20)
fit1 <- coxme(Surv(time, status) ~ age + trt + (trt|inst) + strata(inst),
              vfixed=.3, simdata, refine.n=nsim, refine.detail=TRUE)

debug <- fit1$refine.detail 

nfrail <- length(unlist(ranef(fit1)))
hmatbb <- fit1$hmat[1:nfrail, 1:nfrail]
bhat <- unlist(ranef(fit1))  #random coefs
set.seed(20)
rdf <- coxme.control()$refine.df
bmat <- matrix(rnorm(nfrail*nsim), nfrail)  # replicate the simulations
bmat <- backsolve(hmatbb, bmat) /
    rep(sqrt(rchisq(nsim, rdf)/rdf), each=nfrail)
bmat <- bmat + bhat
if (!is.null(debug)) all.equal(bmat, debug$bmat)

clog <- double(nsim)
Xmat <- scale(as.matrix(simdata[,c('age', 'trt')]), fit1$means, FALSE)

# Part 1, loglik for a set of nearby Cox models
fix.lp <- Xmat %*% fixef(fit1)
for (i in 1:nsim) {
    lp <- fix.lp + bmat[simdata$inst,i]*simdata$trt
    tfit <- coxph(Surv(time, status) ~ offset(lp) + strata(inst), simdata)
    clog[i] <- tfit$loglik
    }
if (!is.null(debug)) aeq(clog, debug$loglik)

# Part 2: Taylor series for the PPL
b.sig <- t(bmat-bhat) %*% hmatbb  #b time sqrt(H)
taylor <- rowSums(b.sig^2)/2

temp2 <- cbind(clog-colSums(bmat^2)/.6 , fit1$log[3] - taylor)
m2 <- fit1$log[2]
errhat <- exp(temp2[,1]-m2) - exp(temp2[,2]-m2)

require(mvtnorm)
tdens <- dmvt(t(bmat), delta=bhat, sigma=as.matrix(solve(hmatbb)), df=rdf)
gdens <- -(log(2*pi) + log(.3))*nfrail/2

errhat <- errhat * exp(gdens-tdens)  # account for importance sampling dis
if (!is.null(debug)) {
    aeq(errhat, debug$errhat)
    }
mtemp <- mean(errhat)
aeq(c(log(1+ mean(errhat)), sqrt(var(errhat)/nsim)/(1+mtemp)), fit1$refine)

