### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {m,mgf}norm functions to compute raw and the
### moment generating function for the Normal distribution (as defined
### in R).
###
### See Chapter 13 of Johnson & Kotz, Continuous univariate
### distributions, volume 1, Wiley, 1970
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mnorm <- function(order, mean = 0, sd = 1)
    .External("actuar_do_dpq", "mnorm", order, mean, sd, FALSE)

mgfnorm <- function(x, mean = 0, sd = 1, log = FALSE)
    .External("actuar_do_dpq", "mgfnorm", x, mean, sd, log)
