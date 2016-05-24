#############################################################################
##                                                                         ##
##   Tests for methods                                                     ##
##                                                                         ##
#############################################################################

## --- Load test routines and test parameters -------------------------------

source("test_routines.R")

## --- CONT: Chi^2 goodness-of-fit test -------------------------------------

## TDR (Transformed Density Rejection)

adv.tdr.norm <- function (n) {
        pdf <- function (x) { exp(-0.5*x^2) }
        dpdf <- function (x) { -x*exp(-0.5*x^2) }
        dist <- new("unuran.cont", pdf=pdf, dpdf=dpdf, lb=-Inf, ub=Inf)
        gen <- unuran.new(dist, "tdr")
        ur(gen,n)
}
unur.test.cont("adv.tdr.norm", rfunc=adv.tdr.norm, pfunc=pnorm)
rm(adv.tdr.norm)

adv.tdr.norm.wl <- function (n) {
        logpdf <- function (x) { -0.5*x^2 }
        dlogpdf <- function (x) { -x }
        dist <- new("unuran.cont", pdf=logpdf, dpdf=dlogpdf, islog=TRUE, lb=-Inf, ub=Inf)
        gen <- unuran.new(dist, "tdr")
        ur(gen,n)
}
unur.test.cont("adv.tdr.norm.wl", rfunc=adv.tdr.norm.wl, pfunc=pnorm)
rm(adv.tdr.norm.wl)

## --- DISCR: Chi^2 goodness-of-fit test ------------------------------------

## DGT (Discrete Guide Table method)
size <- 100
prob <- 0.3
binom.pmf <- function (x) { dbinom(x, size, prob) }
binom.probs <- dbinom(0:size, size, prob)
adv.dgt.binom <- function (n,lb=0,ub=size) {
        dist <- new("unuran.discr", pv=binom.probs, lb=lb)
        gen <- unuran.new(dist, "dgt")
        ur(gen,n)
}
unur.test.discr("adv.dgt.binom", rfunc=adv.dgt.binom, dfunc=binom.pmf, domain=c(0,size))
unur.test.discr("adv.dgt.binom", rfunc=adv.dgt.binom, pv=binom.probs, domain=c(0,size))
rm(adv.dgt.binom)
rm(size,prob,binom.pmf,binom.probs)

## DAU (Discrete Alias-Urn method)
size <- 100
prob <- 0.3
binom.pmf <- function (x) { dbinom(x, size, prob) }
binom.probs <- dbinom(0:size, size, prob)
adv.dau.binom <- function (n,lb=0,ub=size) {
        dist <- new("unuran.discr", pv=binom.probs, lb=lb)
        gen <- unuran.new(dist, "dau")
        ur(gen,n)
}
unur.test.discr("adv.dau.binom", rfunc=adv.dau.binom, dfunc=binom.pmf, domain=c(0,size))
unur.test.discr("adv.dau.binom", rfunc=adv.dau.binom, pv=binom.probs, domain=c(0,size))
rm(adv.dau.binom)
rm(size,prob,binom.pmf,binom.probs)


## --- CMV: Chi^2 goodness-of-fit test --------------------------------------

samplesize <- 1.e4

## HITRO (Hit-and-Run + Ratio-of-Uniforms)
adv.hitro.norm <- function (n) {
        pdf <- function (x) { exp(-0.5*sum(x^2)) }
        dist <- new("unuran.cmv", dim=2, pdf=pdf)
        gen <- unuran.new(dist, "hitro; thinning=10")
        ur(gen,n)
}
unur.test.cmv("hitro.norm", rfunc=adv.hitro.norm, pfunc=pnorm)
rm(adv.hitro.norm)

## VNROU (Naive  multivariate Ratio-of-Uniforms method)
adv.vnrou.norm <- function (n) {
        pdf <- function (x) { exp(-0.5*sum(x^2)) }
        dist <- new("unuran.cmv", dim=2, pdf=pdf)
        gen <- unuran.new(dist, "vnrou")
        ur(gen,n)
}
unur.test.cmv("vnrou.norm", rfunc=adv.vnrou.norm, pfunc=pnorm)
rm(adv.vnrou.norm)


## -- Print statistics ------------------------------------------------------

unur.test.statistic()

## -- End -------------------------------------------------------------------

detach("package:Runuran",unload = TRUE)

## --------------------------------------------------------------------------
