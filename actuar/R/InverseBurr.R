### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}burr functions to compute
### characteristics of the Inverse Burr distribution. The version used
### in these functions has cumulative distribution function
###
###   Pr[X <= x] = (u/(1 + u))^shape1, u = (x/scale)^shape2, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvburr <- function (x, shape1, shape2, rate = 1, scale = 1/rate,
                      log = FALSE)
    .External("actuar_do_dpq", "dinvburr", x, shape1, shape2, scale, log)

pinvburr <- function(q, shape1, shape2, rate = 1, scale = 1/rate,
                     lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pinvburr", q, shape1, shape2, scale,
              lower.tail, log.p)

qinvburr <- function(p, shape1, shape2, rate = 1, scale = 1/rate,
                     lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qinvburr", p, shape1, shape2, scale,
              lower.tail, log.p)

rinvburr <- function(n, shape1, shape2, rate = 1, scale = 1/rate)
    .External("actuar_do_random", "rinvburr", n, shape1, shape2, scale)

minvburr <- function(order, shape1, shape2, rate = 1, scale = 1/rate)
    .External("actuar_do_dpq", "minvburr", order, shape1, shape2, scale, FALSE)

levinvburr <- function(limit, shape1, shape2, rate = 1, scale = 1/rate,
                       order = 1)
    .External("actuar_do_dpq", "levinvburr", limit, shape1, shape2, scale,
              order, FALSE)
