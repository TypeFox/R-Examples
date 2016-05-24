### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}invtrgamma functions to compute
### characteristics of the Inverse Transformed Gamma distribution. The
### version used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - pgamma((x/scale)^shape2, shape1, scale = 1)
###
### or, equivalently,
###
###   Pr[X <= x] = 1 - pgamma(x^shape2, shape1, scale = scale^shape2).
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvtrgamma <- function (x, shape1, shape2, rate = 1, scale = 1/rate,
                         log = FALSE)
    .External("actuar_do_dpq", "dinvtrgamma", x, shape1, shape2, scale, log)

pinvtrgamma <- function(q, shape1, shape2, rate = 1, scale = 1/rate,
                        lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pinvtrgamma", q, shape1, shape2, scale,
              lower.tail, log.p)

qinvtrgamma <- function(p, shape1, shape2, rate = 1, scale = 1/rate,
                        lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qinvtrgamma", p, shape1, shape2, scale,
              lower.tail, log.p)

rinvtrgamma <- function(n, shape1, shape2, rate = 1, scale = 1/rate)
    .External("actuar_do_random", "rinvtrgamma", n, shape1, shape2, scale)

minvtrgamma <- function(order, shape1, shape2, rate = 1, scale = 1/rate)
    .External("actuar_do_dpq", "minvtrgamma", order, shape1, shape2, scale, FALSE)

levinvtrgamma <- function(limit, shape1, shape2, rate = 1, scale = 1/rate,
                          order = 1)
    .External("actuar_do_dpq", "levinvtrgamma", limit, shape1, shape2, scale,
              order, FALSE)
