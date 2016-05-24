library(coxme)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
#
# Ensure that a call with two terms and a call with one term but zero
#  correlation give the same answer.  These are identical mathmatically,
#  but follow two different indexing paths in the code.
#
approx <- "efron"
set.seed(3101953)
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


# The true model has per site coefficients of
#    site 1:  .5 + .25*trt1
#    site 2: 1.5 + .35*trt1
#    site 3: 2.0 + .45*trt1
#    site 4: 1.0 + .55*trt1
smdata <- mkdata(150)  # 150 per site

# Variance function that sets the parameters to a fixed value
myvar <- function(var=c(.25, .025, .078)) {
    initialize <- function(vinit, fixed, intercept, G, X, sparse, ...) {
        imap <- as.matrix(as.numeric(G[[1]]))
        list(theta=NULL, imap=imap, X=X, xmap=imap+4,
             parms=list(theta=var, fixed=rep(TRUE,8),
                        xname=dimnames(X)[[2]], nvar=ncol(X)))
    }
    generate <- function(newtheta, parms) {
        theta <- parms$theta
        varmat <- diag(rep(theta[1:2], each=4))
        for (i in 1:4) varmat[i,i+4] <- varmat[i+4,i] <- theta[3]
        varmat
    }
    
    wrapup <- function(newtheta, b, parms) {
        theta <- parms$theta
        rtheta <- list(site=c(var1=theta[1], var2=theta[2], covar=theta[3]))
        b <- matrix(b, ncol=2)
        dimnames(b) <- list(1:4, c("Intercept", "Slope"))
        list(theta=rtheta, b=b)
    }
    out <- list(initialize=initialize, generate=generate, wrapup=wrapup)
    class(out) <- 'coxmevar'
    out
    }

fit1<-  coxme(Surv(futime, status) ~ trt1 + (1 | site) + (trt1|site), smdata,
              ties=approx)
fit2 <- coxme(Surv(futime, status) ~ trt1 + (1+trt1 | site), smdata,
              ties=approx, varlist=myvar(c(unlist(VarCorr(fit1)), 0)))

aeq(fit1$log, fit2$log)
aeq(unlist(ranef(fit1)), unlist(ranef(fit2)))
aeq(fixef(fit1), fixef(fit2))

# Same models, but fit with the start,stop part of the code
dummy <- runif(nrow(smdata), -4, -1)  #all start times before any stop times
fit1b <- coxme(Surv(dummy, futime, status) ~ trt1 + (1 | site) + (trt1|site), 
               smdata, ties=approx)
all.equal(ranef(fit1b), ranef(fit1))
aeq(fit1$loglik, fit1b$loglik)

fit2b <- coxme(Surv(dummy, futime, status) ~ trt1 + (1+trt1 | site), smdata,
              ties=approx, varlist=myvar(c(unlist(VarCorr(fit1)), 0)))
all.equal(fit2b$loglik, fit2$loglik)
all.equal(coef(fit2b), coef(fit2))
