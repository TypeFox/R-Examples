### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {m,lev,mgf}chisq functions to compute raw and
### limited moments, and the moment generating function for
### the Chi-square distribution (as defined in R)
###
### See Chapter 17 of Johnson & Kotz, Continuous univariate
### distributions, volume 1, Wiley, 1970
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mchisq <- function(order, df, ncp = 0)
    .External("actuar_do_dpq", "mchisq", order, df, ncp, FALSE)

levchisq <- function(limit, df, ncp = 0, order = 1)
    .External("actuar_do_dpq", "levchisq", limit, df, ncp, order, FALSE)

mgfchisq <- function(x, df, ncp = 0, log = FALSE)
    .External("actuar_do_dpq", "mgfchisq", x, df, ncp, log)
