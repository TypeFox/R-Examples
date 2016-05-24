### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}invparalogis functions to compute
### characteristics of the Inverse Paralogistic distribution. The
### version used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = (u/(1 + u))^shape, u = (x/scale)^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvparalogis <- function (x, shape, rate = 1, scale = 1/rate, log = FALSE)
    .External("actuar_do_dpq", "dinvparalogis", x, shape, scale, log)

pinvparalogis <- function(q, shape, rate = 1, scale = 1/rate,
                          lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pinvparalogis", q, shape, scale, lower.tail, log.p)

qinvparalogis <- function(p, shape, rate = 1, scale = 1/rate,
                          lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qinvparalogis", p, shape, scale, lower.tail, log.p)

rinvparalogis <- function(n, shape, rate = 1, scale = 1/rate)
    .External("actuar_do_random", "rinvparalogis", n, shape, scale)

minvparalogis <- function(order, shape, rate = 1, scale = 1/rate)
     .External("actuar_do_dpq", "minvparalogis", order, shape, scale, FALSE)

levinvparalogis <- function(limit, shape, rate = 1, scale = 1/rate,
                            order = 1)
     .External("actuar_do_dpq", "levinvparalogis", limit, shape, scale, order, FALSE)
