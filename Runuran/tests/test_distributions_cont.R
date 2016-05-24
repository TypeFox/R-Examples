#############################################################################
##                                                                         ##
##   Tests for wrapper functions for special distributions                 ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Continuous univariate distributions                                   ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Remark: You must use named arguments when calling the test routines!  ##
##                                                                         ##
#############################################################################

## --- Load test routines and test parameters -------------------------------

source("test_routines.R")

## --- Chi^2 goodness-of-fit test -------------------------------------------

## Beta distribution - (replacement for rbeta) ------------------------------
for (i in 1:n.rep.domains)
        unur.test.cont("beta", shape1=runif(1,0.2,10), shape2=runif(1,0.2,10), 
                       domain=sort(runif(2)))
for (i in 1:n.rep.params)
        unur.test.cont("beta", shape1=runif(1,0.2,10), shape2=runif(1,0.2,10))

## Burr family of distributions ---------------------------------------------

## Cauchy distribution - (replacement for rcauchy) --------------------------
unur.test.cont("cauchy")
for (i in 1:n.rep.domains)
        unur.test.cont("cauchy", domain=sort(rnorm(2)))
for (i in 1:n.rep.params)
        unur.test.cont("cauchy", location=rcauchy(1), scale=rgamma(1,shape=2))

## Chi distribution ---------------------------------------------------------

## Chi^2 distribution - (replacement for rchisq) ----------------------------
for (i in 1:n.rep.domains)
        unur.test.cont("chisq", df=3, domain=sort(rexp(2)))
for (i in 1:n.rep.params)
        unur.test.cont("chisq", df=rgamma(1,shape=2))

## Exponential distribution - (replacement for rexp) ------------------------
for (i in 1:n.rep.domains)
        unur.test.cont("exp", domain=sort(rexp(2)))
for (i in 1:n.rep.params)
        unur.test.cont("exp", rate=runif(1,0.1,10))

## Extreme value type I (Gumbel-type) distribution --------------------------

## Extreme value type II (Frechet-type) distribution ------------------------

## F distribution  - (replacement for rf) -----------------------------------
for (i in 1:n.rep.domains)
        unur.test.cont("f", df1=runif(1,0.1,10), df2=runif(1,0.1,10), 
                       domain=sort(runif(2)))
for (i in 1:n.rep.params)
        unur.test.cont("f", df1=runif(1,0.1,10), df2=runif(1,0.1,10))

## Gamma distribution  - (replacement for rgamma) ---------------------------
for (i in 1:n.rep.domains)
        unur.test.cont("gamma", shape=runif(1,0.1,10), domain=sort(rexp(2)))
for (i in 1:n.rep.params)
        unur.test.cont("gamma", shape=runif(1,0.1,10), scale=runif(1,0.1,10))

## Generalized inverse Gaussian ---------------------------------------------

## Hyperbolic distribution --------------------------------------------------

## Laplace (double exponential) distribution --------------------------------

## Lognormal distribution  - (replacement for rlnorm) -----------------------
unur.test.cont("lnorm")
for (i in 1:n.rep.domains)
        unur.test.cont("lnorm", domain=sort(rexp(2)))
for (i in 1:n.rep.params)
        unur.test.cont("lnorm", meanlog=rnorm(1), sdlog=rgamma(1,shape=2))

## Logistic distribution - (replacement for rlogistic) ----------------------
unur.test.cont("logis")
for (i in 1:n.rep.domains)
        unur.test.cont("logis", domain=sort(rnorm(2)))
for (i in 1:n.rep.params)
        unur.test.cont("logis", location=rcauchy(1), scale=rgamma(1,shape=2))

## Lomax distribution (Pareto distribution of second kind) ------------------

## Normal (Gaussian) distribution - (replacement for rnorm) -----------------
unur.test.cont("norm")
for (i in 1:n.rep.domains)
        unur.test.cont("norm", domain=sort(rnorm(2)))
for (i in 1:n.rep.params)
        unur.test.cont("norm", mean=rcauchy(1), sd=rgamma(1,shape=2))

rud <- function (n) {
   dist <- udnorm()
   gen <- unuran.new(dist)
   ur(gen,n)
}
unur.test.cont("udnorm", rfunc=rud, pfunc=pnorm)
rm(rud)

## Pareto distribution ------------------------------------------------------

## Planck distribution ------------------------------------------------------

## Powerexponential (Subbotin) distribution ---------------------------------

## Rayleigh distribution ----------------------------------------------------

## Student's t distribution - (replacement for rt) --------------------------
for (i in 1:n.rep.domains)
        unur.test.cont("t", df=2, domain=sort(rnorm(2)))
for (i in 1:n.rep.params)
        unur.test.cont("t", df=0.1+rgamma(1,shape=2))

## Weibull distribution - (replacement for rweibull) ------------------------
for (i in 1:n.rep.domains) {
        s <- runif(1,0.1,10)
        unur.test.cont("weibull", shape=s, domain=sort(urweibull(n=2,shape=s)))
}
for (i in 1:n.rep.params)
        unur.test.cont("weibull", shape=runif(1,0.1,10), scale=runif(1,0.1,10))


## -- Print statistics ------------------------------------------------------

unur.test.statistic()

## -- End -------------------------------------------------------------------

detach("package:Runuran",unload = TRUE)

## --------------------------------------------------------------------------
