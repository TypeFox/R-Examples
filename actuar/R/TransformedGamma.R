### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}trgamma functions to compute
### characteristics of the Transformed Gamma distribution. The version
### used in these functions has cumulative distribution function
###
###   Pr[X <= x] = pgamma((x/scale)^shape2, shape1, scale = 1), x > 0
###
### or, equivalently,
###
###   Pr[X <= x] = pgamma(x^shape2, shape1, scale = scale^shape2), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dtrgamma <- function (x, shape1, shape2, rate = 1, scale = 1/rate,
                      log = FALSE)
    .External("actuar_do_dpq", "dtrgamma", x, shape1, shape2, scale, log)

ptrgamma <- function(q, shape1, shape2, rate = 1, scale = 1/rate,
                     lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "ptrgamma", q, shape1, shape2, scale,
              lower.tail, log.p)

qtrgamma <- function(p, shape1, shape2, rate = 1, scale = 1/rate,
                     lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qtrgamma", p, shape1, shape2, scale,
              lower.tail, log.p)

rtrgamma <- function(n, shape1, shape2, rate = 1, scale = 1/rate)
    .External("actuar_do_random", "rtrgamma", n, shape1, shape2, scale)

mtrgamma <- function(order, shape1, shape2, rate = 1, scale = 1/rate)
    .External("actuar_do_dpq", "mtrgamma", order, shape1, shape2, scale, FALSE)

levtrgamma <- function(limit, shape1, shape2, rate = 1, scale = 1/rate,
                       order = 1)
    .External("actuar_do_dpq", "levtrgamma", limit, shape1, shape2, scale, order, FALSE)
