### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}paralogis functions to compute
### characteristics of the paralogistic distribution. The version used
### in these functions has cumulative distribution function
###
###   Pr[X <= x] = 1 - (1/(1 + (x/scale)^shape))^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dparalogis <- function (x, shape, rate = 1, scale = 1/rate, log = FALSE)
    .External("actuar_do_dpq", "dparalogis", x, shape, scale, log)

pparalogis <- function(q, shape, rate = 1, scale = 1/rate,
                       lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pparalogis", q, shape, scale, lower.tail, log.p)

qparalogis <- function(p, shape, rate = 1, scale = 1/rate,
                       lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qparalogis", p, shape, scale, lower.tail, log.p)

rparalogis <- function(n, shape, rate = 1, scale = 1/rate)
    .External("actuar_do_random", "rparalogis", n, shape, scale)

mparalogis <- function(order, shape, rate = 1, scale = 1/rate)
    .External("actuar_do_dpq", "mparalogis", order, shape, scale, FALSE)

levparalogis <- function(limit, shape, rate = 1, scale = 1/rate,
                         order = 1)
    .External("actuar_do_dpq", "levparalogis", limit, shape, scale, order, FALSE)
