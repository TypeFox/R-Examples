#### Tests for "Functionals":  Root finding, Optimization, Integration, etc

stopifnot(require("Rmpfr"))

(f.chk <- system.file("check-tools.R", package="Rmpfr", mustWork=TRUE))
source(f.chk, keep.source=FALSE)
## -> all.eq.finite(), all.EQ(), ...

options(warn = 1)# warnings *immediately*
(doExtras <- Rmpfr:::doExtras())

### 1. Integration -----------------------------------------------

## Example from Lauren K, June 2014 (~/R/MM/Pkg-ex/Rmpfr/integrateR-LaurenK.R):
beta0  <- 0.05
beta1  <- 0.05
(tau <- sqrt(0.01*0.05*0.95/0.99))# = 0.0219..
##
Z00  <- 9
Z01  <- 1
Z10  <- 18
Z11  <- 2
N <- Z00+Z01+Z10+Z11

integrand <- function(u) {
    ee.u <- exp(-u^2/2)/(sqrt(2*pi)*tau)
    b0u <- beta0 + tau*u
    b1u <- beta1 + b0u # == beta0+beta1+ tau*u

    ee.u ^ (Z00+Z01 + Z10+Z11) *
        (1-b0u)^Z00 * b0u ^ Z01 *
        (1-b1u)^Z10 * b1u ^ Z11
}

## MM: note how the integrand() function looks:
op <- par(mfcol=c(3,1), mgp=c(1.5, .6, 0), mar=.1+c(3,3,0,0))
  ax1 <- function(a,b) axis(1, at=c(a,b), col=2, col.axis=2, tcl= +3/4, mgp=c(3,0,0))
  curve(integrand, -5,5, n=1000)
  cc <- adjustcolor(2, 1/4) ## look closer:
  ep <- .01; rect(-3, -ep, 0, +ep, col=cc, border=cc); ax1(-3,0)
  curve(integrand, -3,0, n=1000, ylim = c(-ep,ep))
  ## but now look really closely :
  ep <- .001; rect(-3, -ep, -2, +ep, col=cc); ax1(-3,-2)
  curve(integrand, -3,-2, n=1000, ylim = c(-ep, ep))
par(op)

(I1  <- integrate(integrand,lower = -100, upper = 100))
(I1. <- integrate(integrand,lower = -100, upper = 100, rel.tol = 1e-14))

showSys.time(I2 <- integrateR(integrand, lower = -100, upper = 100))
I2 ## ... Warning ‘no convergence up to order  13’
## Solaris Sparc (2014-06, CRAN checks); thanks Brian: print(I2[1:2], digits=15)
I2.Solaris <- list(value = 1.3963550396006e+33, abs.error = 1.79487857486724e+28)
I.db <- list(value = 1.39635503960059e+33, abs.error = 1.79487857478077e+28)
stopifnot(
    all.equal(I2[1:2], I.db, tol = 1e-10)# Solaris SPARC needs at least 4.8e-11
)

## Now using high accuracy
showSys.time(I3 <- integrateR(integrand, lower = mpfr(-100, precBits=256), upper = 100))
## much slower but not better (and not worse)
I3
assert.EQ.(sapply(I3[1:2], asNumeric), unlist(I.db))
## Really get better when decreasing the integration interval
## from  [-100, 100] to  [-10, 10] ... which should give "the same"
showSys.time(I4 <- integrateR(integrand, lower = mpfr(-10, precBits=256), upper = 10,
                              ord = 15, verbose=TRUE))
## ~ 6.6 sec [lynne 2013]
I4

## on the left side, there is "nothing" (and negative, as we know!):
showSys.time(I0.1 <- integrateR(integrand, lower = mpfr(-1000, precBits=256),
                                upper = -10,  ord= 11, verbose=TRUE))
showSys.time(I0.2 <- integrateR(integrand, lower = mpfr(10, precBits=256),
                                upper = 1000, ord= 11, verbose=TRUE))
I0.1
I0.2
I4

## Integral [-1000, +1000 ] = Int[-1000, -10] + Int[-10, +10] + Int[+10, +1000]:

I4 $value + I0.1 $value + I0.2 $value
## but this is really the same as just the middle:
stopifnot(I4 $value + I0.1 $value + I0.2 $value
          == I4 $value)

value <- I4$value; delta <- I4$abs.err
nDig <- -asNumeric(log10(delta/value))
cat("Correct number of digits: ", round(nDig, 2),"\n",
    "Integral I =              ", format(I4$value, digits = ceiling(nDig)),
    " (last change change= ", format(delta, digits = 7),")\n",
    "integrate(.) =            ", format(I1 $value, digits = 22),"\n",
    "integrate(., rtol=1e-15)= ", format(I1.$value, digits = 22),"\n", sep="")




### 2. Root Finding ----------------------------------------------

### 3. Optimization / Minimization, .. ---------------------------



cat('Time elapsed: ', proc.time(),'\n') # "stats"
