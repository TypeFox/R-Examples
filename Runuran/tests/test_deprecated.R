#############################################################################
##                                                                         ##
##   Tests for deprecated function                                         ##
##                                                                         ##
#############################################################################

## --- Load test routines and test parameters -------------------------------

source("test_routines.R")

## --- DISCR: Chi^2 goodness-of-fit test ------------------------------------

## DGT (Discrete Guide Table method)
size <- 100
prob <- 0.3
binom.pmf <- function (x) { dbinom(x, size, prob) }
binom.probs <- dbinom(0:size, size, prob)
urdgt.binom <- function (n,lb=0,ub=size) {
        urdgt(n, probvector=binom.probs)
}
unur.test.discr("urdgt.binom", rfunc=urdgt.binom, dfunc=binom.pmf, domain=c(0,size))
unur.test.discr("urdgt.binom", rfunc=urdgt.binom, pv=binom.probs, domain=c(0,size))
rm(urdgt.binom)
rm(size,prob,binom.pmf,binom.probs)

## DAU (Discrete Alias-Urn method)
size <- 100
prob <- 0.3
binom.pmf <- function (x) { dbinom(x, size, prob) }
binom.probs <- dbinom(0:size, size, prob)
urdau.binom <- function (n,lb=0,ub=size) {
        urdau(n, probvector=binom.probs)
}
unur.test.discr("urdau.binom", rfunc=urdau.binom, dfunc=binom.pmf, domain=c(0,size))
unur.test.discr("urdau.binom", rfunc=urdau.binom, pv=binom.probs, domain=c(0,size))
rm(urdau.binom)
rm(size,prob,binom.pmf,binom.probs)

## -- Print statistics ------------------------------------------------------

unur.test.statistic()

## -- End -------------------------------------------------------------------

detach("package:Runuran",unload = TRUE)

## --------------------------------------------------------------------------
