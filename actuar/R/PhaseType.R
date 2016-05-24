### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,r,m,mgf}ph functions to compute
### characteristics of Phase-type distributions with cumulative
### distribution function
###
###       Pr[X <= x] = 1 - pi %*% exp(Tx) %*% e,
###
### where 'pi' is the initial probability vector, 'T' is the
### subintensity matrix and 'e' is 1-vector of R^m.
###
### See Bladt, M. (2005), "A review on phase-type distributions and
### their use in risk theory", Astin Bulletin 35, p. 145-161.
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dphtype <- function (x, prob, rates, log = FALSE)
    .External("actuar_do_dpqphtype", "dphtype", x, prob, rates, log)

pphtype <- function (q, prob, rates, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpqphtype", "pphtype", q, prob, rates, lower.tail, log.p)

rphtype <- function (n, prob, rates)
    .External("actuar_do_randomphtype", "rphtype", n, prob, rates)

mphtype <- function (order, prob, rates)
    .External("actuar_do_dpqphtype", "mphtype", order, prob, rates, FALSE)

mgfphtype <- function (x, prob, rates, log = FALSE)
    .External("actuar_do_dpqphtype", "mgfphtype", x, prob, rates, log)
