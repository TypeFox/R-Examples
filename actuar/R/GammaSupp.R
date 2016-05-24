### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {m,lev,mgf}gamma functions to compute raw and
### limited moments, and the moment generating function for
### the Gamma distribution (as defined in R)
###
### See Chapter 17 of Johnson & Kotz, Continuous univariate
### distributions, volume 1, Wiley, 1970
###
### AUTHORS: Mathieu Pigeon, Christophe Dutang,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

mgamma <- function(order, shape, rate = 1, scale = 1/rate)
    .External("actuar_do_dpq", "mgamma", order, shape, scale, FALSE)

levgamma <- function(limit, shape, rate = 1, scale = 1/rate, order = 1)
    .External("actuar_do_dpq", "levgamma", limit, shape, scale, order, FALSE)

mgfgamma <- function(x, shape, rate = 1, scale = 1/rate, log = FALSE)
    .External("actuar_do_dpq", "mgfgamma", x, shape, scale, log)
