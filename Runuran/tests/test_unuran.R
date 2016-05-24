
#############################################################################
##                                                                         ##
##   Tests for class 'unuran'                                              ##
##                                                                         ##
#############################################################################

## --- Test Parameters ------------------------------------------------------

## size of sample for test
samplesize <- 1.e6

## break points for chi^2 GoF test
nbins <- as.integer(sqrt(samplesize))
breaks <- (0:nbins)/nbins

## level of significance
alpha <- 1.e-3

## seed for uniform RNG
set.seed(123456)

## --- Load library ---------------------------------------------------------

library(Runuran)


#############################################################################
##                                                                          #
##  Auxiliary routines                                                      #
##                                                                          #
#############################################################################

## -- Auxiliary routines ------------------------------------------------------

## Test whether there is an error -------------------------------------------
is.error <- function (expr) { is(try(expr), "try-error") }


## --- Continuous distributions ---------------------------------------------

## Create an object
unr <- new("unuran", "normal()")

## Print object
unr
print(unr)
unuran.details(unr)
unuran.details(unr,show=TRUE,return.list=FALSE)
unuran.details(unr,show=TRUE,return.list=TRUE)
unuran.details(unr,show=FALSE,return.list=TRUE)
print(unuran.details(unr,show=FALSE,return.list=TRUE))
unuran.details(unr,show=FALSE,return.list=FALSE)
print(unuran.details(unr,show=FALSE,return.list=FALSE))

## Test object properties
unuran.is.inversion(unr)

## Draw samples
unuran.sample(unr)
unuran.sample(unr,10)
x <- unuran.sample(unr, samplesize)
ur(unr)
ur(unr,10)

## Run a chi-square GoF test
chisq.test( hist(pnorm(x),plot=FALSE,breaks=breaks)$density )


## Create an object
unr <- unuran.new("normal()")

## Draw samples
unuran.sample(unr)
unuran.sample(unr,10)
x <- unuran.sample(unr, samplesize)

## Run a chi-square GoF test
pval <- chisq.test( hist(pnorm(x),plot=FALSE,breaks=breaks)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))

## another example (for testing print)
unr <- new("unuran", "normal()", "pinv")
unr
print(unr)
unuran.details(unr)
print(unuran.details(unr,show=TRUE,return.list=TRUE))
unuran.is.inversion(unr)

## remove (so that valgrind does not see lost memory from UNU.RAN)
rm(unr)


## --- Continuous distributions - S4 distribution object --------------------

## use PDF
gausspdf <- function (x) { exp(-0.5*x^2) }
gaussdpdf <- function (x) { -x*exp(-0.5*x^2) }
gauss <- new("unuran.cont", pdf=gausspdf, dpdf=gaussdpdf, lb=-Inf, ub=Inf, center=0.1)
unr <- 0; unr <- unuran.new(gauss, "tdr")
unr
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pnorm(x),plot=FALSE,breaks=breaks)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))
rm(unr)

## use PDF
gausspdf <- function (x) { exp(-0.5*x^2) }
gaussdpdf <- function (x) { -x*exp(-0.5*x^2) }
gauss <- unuran.cont.new(pdf=gausspdf, dpdf=gaussdpdf, lb=-Inf, ub=Inf, center=0.1)
unr <- unuran.new(gauss, "tdr")
unr
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pnorm(x),plot=FALSE,breaks=breaks)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))
rm(unr)

## use logPDF
gausspdf <- function (x) { -0.5*x^2 }
gaussdpdf <- function (x) { -x }
gauss <- new("unuran.cont", pdf=gausspdf, dpdf=gaussdpdf, islog=TRUE, lb=-Inf, ub=Inf, mode=0)
unr <- unuran.new(gauss, "tdr")
unr
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pnorm(x),plot=FALSE,breaks=breaks)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))
rm(unr)

## use logPDF (use ARS to test 'print' function)
gausspdf <- function (x) { -0.5*x^2 }
gaussdpdf <- function (x) { -x }
gauss <- new("unuran.cont", pdf=gausspdf, dpdf=gaussdpdf, islog=TRUE, lb=-Inf, ub=Inf)
unr <- unuran.new(gauss, "ars")
unr
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pnorm(x),plot=FALSE,breaks=breaks)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))
rm(unr)


## --- Discrete distributions -----------------------------------------------

## Create an object
unr <- new("unuran", "binomial(20,0.5)", "dgt")

## Draw samples
unuran.sample(unr)
unuran.sample(unr,10)
x <- unuran.sample(unr, samplesize)
rm(unr)


## --- Discrete distributions - S4 distribution object ----------------------

## use PV
pv <- dbinom(0:100,100,0.3)
binom <- new("unuran.discr",pv=pv,lb=0)
unr <- unuran.new(binom, "dgt")
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pbinom(x,100,0.3),plot=FALSE)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))
rm(unr)

## use PMF
pmf <- function(x) dbinom(x,100,0.3)
binom <- new("unuran.discr",pmf=pmf,lb=0,ub=100)
unr <- unuran.new(binom, "dgt")
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pbinom(x,100,0.3),plot=FALSE)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))
rm(unr)

## use PMF
pmf <- function(x) dbinom(x,100,0.3)
binom <- new("unuran.discr",pmf=pmf,lb=0,ub=100)
unr <- unuran.new(binom, "dari")
x <- unuran.sample(unr, samplesize)
pval <- chisq.test( hist(pbinom(x,100,0.3),plot=FALSE)$density )$p.value
if (pval < alpha) stop("chisq test FAILED!  p-value=",signif(pval))
rm(unr)


## --- Continuous Multivariate distributions --------------------------------

mvpdf <- function (x) { exp(-sum(x^2)) }
mvd <- new("unuran.cmv", dim=2, pdf=mvpdf, mode=c(0,0))
unr <- unuran.new(mvd, "hitro")
x <- unuran.sample(unr, 10)
x
rm(unr)

unr <- unuran.new(mvd, "vnrou")
x <- unuran.sample(unr, 10)
x
rm(unr)


## --- quantile function ----------------------------------------------------

## test U-error
unr <- unuran.new("normal()","hinv; u_resolution=1.e-13")

## single double as argument
Tmax <- 0
for (U in (0:20)/20) {
        T <- pnorm( uq(unr,U) ) - U
        print (T)
        Tmax <- max(abs(T),Tmax)
}
cat("Max. error =",Tmax,"\n")
if (Tmax > 1.e-13) stop ("Max. error exceeds limit")

## vector argument
U <- (0:20)/20
T <- pnorm( uq(unr,U) ) - U
print (T)
Tmax <- max(abs(T))
cat("Max. error =",Tmax,"\n")
if (Tmax > 1.e-13) stop ("Max. error exceeds limit")

## special arguments
for (U in c(-1,-0.001,0,0.5,1.,1.001,NA,NaN) ) {
        cat ("U =",U,"\tX =",uq(unr,U),"\n")
}

U <- c(-1,-0.001,0,0.5,1.,1.001,NA,NaN)
T <- uq(unr,U)
rbind(U,T)

uq(unr,numeric())

rm(unr)

## test whether 'uq' throws an error when UNU.RAN object does not implement
## an inversion method
unr <- unuran.new("normal()","tdr")
if( ! is.error( uq(unr,0.5) ) )
   stop("'uq' does not detect UNU.RAN object with non-inversion method")

## test whether 'uq' detects invalid arguments
if( ! is.error( uq(1,0.5) ) )
   stop("'uq' does not detect invalid argument 'unr'")
        
if( ! is.error( uq(unr,"a") ) )
   stop("'uq' does not detect invalid argument 'U'")
        
rm(unr)


## --- density function -----------------------------------------------------

distr <- unuran.cont.new(pdf=function(x){exp(-x)}, lb=0,ub=Inf)
x <- rexp(100)
e <- max(abs(ud(distr, x) - dexp(x)))
e; if (e>1.e-10) stop("error too large")

distr <- unuran.cont.new(pdf=function(x){-x}, islog=TRUE, lb=0,ub=Inf)
x <- rexp(100)
e <- max(abs(ud(distr, x, islog=TRUE) + x))
e; if (e>1.e-10) stop("error too large")

distr <- udgeom(prob=0.8)
x <- rgeom(100, prob=0.8)
e <- max(abs(ud(distr,x) - dgeom(x,prob=0.8)))
e; if (e>1.e-10) stop("error too large")

rm(distr,x,e)

unr <- pinv.new(pdf=function(x){exp(-x)}, lb=0,ub=Inf)
x <- rexp(100)
e <- max(abs(ud(unr, x) - dexp(x)))
e; if (e>1.e-10) stop("error too large")

unr <- pinv.new(pdf=function(x){-x}, islog=TRUE, lb=0,ub=Inf)
x <- rexp(100)
e <- max(abs(ud(unr, x, islog=TRUE) +x))
e; if (e>1.e-10) stop("error too large")

unr <- darid.new(udgeom(prob=0.8))
x <- rgeom(100, prob=0.8)
e <- max(abs(ud(unr,x) - dgeom(x,prob=0.8)))
e; if (e>1.e-10) stop("error too large")

rm(unr,x,e)


distr <- unuran.cont.new(lb=0,ub=1)
if (!all(is.na(ud(distr,1))))
  stop("'ud' ignores missing PDF")

distr <- unuran.discr.new(lb=0,ub=1)
if (!all(is.na(ud(distr,1))))
  stop("'ud' ignores missing PMF")

unr <- pinv.new(cdf=pexp,lb=0,ub=Inf)
if (!all(is.na(ud(distr,1))))
  stop("'ud' ignores missing PDF")

unr <- pinv.new(pdf=dexp,lb=0,ub=Inf)
unuran.packed(unr) <- TRUE
if( ! is.error( ud(unr, 1) ) ) 
  stop("'ud' does not detect packed generator object")

rm(distr,unr)


## --- distribution function ------------------------------------------------

distr <- unuran.cont.new(cdf=function(x){1-exp(-x)}, lb=0,ub=Inf)
x <- rexp(100)
e <- max(abs(up(distr, x) - pexp(x)))
e; if (e>1.e-10) stop("error too large")

distr <- udgeom(prob=0.8)
x <- rgeom(100, prob=0.8)
e <- max(abs(up(distr,x) - pgeom(x,prob=0.8)))
e; if (e>1.e-10) stop("error too large")

unr <- pinv.new(pdf=function(x){exp(5-x)}, lb=0,ub=Inf)
x <- rexp(100)
e <- max(abs(up(unr, x) - pexp(x)))
e; if (e>1.e-10) stop("error too large")

rm(distr,unr,x,e)


distr <- unuran.cont.new(lb=0,ub=1)
if( ! is.error( up(distr,1) ) )
  stop("'up' ignores missing CDF")

distr <- unuran.discr.new(lb=0,ub=1)
if( ! is.error( up(distr,1) ) )
  stop("'up' ignores missing CDF")

unr <- tdr.new(pdf=dexp,lb=0,ub=Inf)
if( ! is.error( up(unr,1) ) )
  stop("'up' ignores invalid method PINV")

unr <- pinv.new(pdf=dexp,lb=0,ub=Inf)
unuran.packed(unr) <- TRUE
if( ! is.error( up(unr, 1) ) ) 
  stop("'up' does not detect packed generator object")

rm(distr,unr)


## --- pack -----------------------------------------------------------------

## check print with unuran.packed
unr <- pinv.new(dnorm,lb=0,ub=Inf)
unuran.packed(unr) <- TRUE
unr
unuran.details(unr)
unuran.details(unr,show=TRUE,return.list=TRUE)
unuran.details(unr,show=FALSE,return.list=TRUE)
print(unuran.details(unr,show=FALSE,return.list=TRUE))
unuran.details(unr,show=FALSE,return.list=FALSE)
rm(unr)

## check whether un/packing works/fails
unr <- pinv.new(dnorm,lb=0,ub=Inf)
## this should be o.k.
unuran.packed(unr) <- FALSE
## this should fail with error
unuran.packed(unr) <- TRUE
if( ! is.error( unuran.packed(unr) <- FALSE ) )
   stop("'unuran.packed' tries to unpack UNU.RAN objects")
## this should be o.k.
unuran.packed(unr) <- TRUE
rm(unr)

## test whether non-packable objects are treated correctly
unr <- tdr.new(dnorm,lb=0,ub=Inf)
if( ! is.error( unuran.packed(unr) <- TRUE ) )
   stop("'unuran.packed' does not detect non-packable UNU.RAN objects")
rm(unr)


## --- mixture --------------------------------------------------------------

comp <- c(unuran.new("normal"),unuran.new("cauchy"),unuran.new("exponential"))
prob <- c(1,2,3)
unr <- mixt.new(prob,comp)
x <- unuran.sample(unr, 10)
x
rm(unr)

comp <- c(pinvd.new(udnorm(lb=-Inf,ub=-1)),
          pinvd.new(udcauchy(lb=-1,ub=1)),
          pinvd.new(udexp(lb=1,ub=Inf)) )
prob <- c(1,2,3)
unr <- mixt.new(prob,comp,inversion=TRUE)
x <- unuran.sample(unr, 10)
x
rm(unr)


## --- Check generator ------------------------------------------------------

## run with defaults
unr <- tdrd.new(udnorm())
unuran.verify.hat(unr)
rm(unr)

## do not show result
unr <- tdrd.new(udnorm())
unuran.verify.hat(unr,show=FALSE)
rm(unr)

## another example
pdf <- function(x) { -x - 0.999*log(x) }
dpdf <- function(x) { -1 - 0.999/x }
distr <- unuran.cont.new(pdf=pdf, dpdf=dpdf, islog=TRUE, lb=0, ub=0.5, mode=0)
unr <- unuran.new(distr,method="itdr;cp=-0.999")
unuran.verify.hat(unr)
rm(unr); rm(distr); rm(pdf); rm(dpdf)

## example where the hat violates condition hat(x) >= pdf(x)
## it also stores the ratio in variable 'failed'
unr <- unuran.new(udnorm(),method="nrou;v=0.6")
failed <- unuran.verify.hat(unr)
failed
rm(unr); rm (failed)

## Example for a method that does not implement rejection
unr <- pinvd.new(udnorm())
if (! is.error( unuran.verify.hat(unr) ) )
  stop("method must throw error for inversion method")
rm(unr)


## --- End ------------------------------------------------------------------

silent <- gc()
detach("package:Runuran",unload = TRUE)

## --------------------------------------------------------------------------
