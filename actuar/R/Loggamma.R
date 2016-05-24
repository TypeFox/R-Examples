### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r}lgamma functions to compute
### characteristics of the Loggamma distribution. The version used in
### these functions has cumulative distribution function
###
###   Pr[X <= x] = pgamma(log(x), shape = shapelog, rate = ratelog).
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dlgamma <- function(x, shapelog, ratelog, log = FALSE)
    .External("actuar_do_dpq", "dlgamma", x, shapelog, ratelog, log)

plgamma <- function(q, shapelog, ratelog, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "plgamma", q, shapelog, ratelog, lower.tail, log.p)

qlgamma <- function(p, shapelog, ratelog, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qlgamma", p, shapelog, ratelog, lower.tail, log.p)

rlgamma <- function(n, shapelog, ratelog)
    .External("actuar_do_random", "rlgamma", n, shapelog, ratelog)

mlgamma <- function(order, shapelog, ratelog)
    .External("actuar_do_dpq", "mlgamma", order, shapelog, ratelog, FALSE)

levlgamma <- function(limit, shapelog, ratelog, order = 1)
    .External("actuar_do_dpq", "levlgamma", limit, shapelog, ratelog, order, FALSE)
