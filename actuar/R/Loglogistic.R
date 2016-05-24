### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r}llogis functions to compute
### characteristics of the loglogistic distribution. The version used
### in these functions has cumulative distribution function
###
###   Pr[X <= x] = u/(1 + u), u = (x/scale)^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dllogis <- function (x, shape, rate = 1, scale = 1/rate, log = FALSE)
     .External("actuar_do_dpq", "dllogis", x, shape, scale, log)

pllogis <- function(q, shape, rate = 1, scale = 1/rate,
                    lower.tail = TRUE, log.p = FALSE)
     .External("actuar_do_dpq", "pllogis", q, shape, scale, lower.tail, log.p)

qllogis <- function(p, shape, rate = 1, scale = 1/rate,
                    lower.tail = TRUE, log.p = FALSE)
     .External("actuar_do_dpq", "qllogis", p, shape, scale, lower.tail, log.p)

rllogis <- function(n, shape, rate = 1, scale = 1/rate)
     .External("actuar_do_random", "rllogis", n, shape, scale)

mllogis <- function(order, shape, rate = 1, scale = 1/rate)
     .External("actuar_do_dpq", "mllogis", order, shape, scale, FALSE)

levllogis <- function(limit, shape, rate = 1, scale = 1/rate,
                      order = 1)
     .External("actuar_do_dpq", "levllogis", limit, shape, scale, order, FALSE)
