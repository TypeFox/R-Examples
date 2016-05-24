### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}single-parameter pareto functions. The single-parameter
### Pareto distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - (min/x)^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dpareto1 <- function (x, shape, min, log = FALSE)
    .External("actuar_do_dpq", "dpareto1", x, shape, min, log)

ppareto1 <- function(q, shape, min, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "ppareto1", q, shape, min, lower.tail, log.p)

qpareto1 <- function(p, shape, min, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qpareto1", p, shape, min, lower.tail, log.p)

rpareto1 <- function(n, shape, min)
    .External("actuar_do_random", "rpareto1", n, shape, min)

mpareto1 <- function(order, shape, min)
     .External("actuar_do_dpq", "mpareto1", order, shape, min, FALSE)

levpareto1 <- function(limit, shape, min, order = 1)
     .External("actuar_do_dpq", "levpareto1", limit, shape, min, order, FALSE)
