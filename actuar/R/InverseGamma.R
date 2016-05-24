### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev,mgf}invgamma functions to compute
### characteristics of the Inverse Gamma distribution. The version
### used in these functions has cumulative distribution function
###
###   Pr[X <= x] = 1 - pgamma(scale/x, shape, scale = 1)
###
### or, equivalently,
###
###   Pr[X <= x] = 1 - pgamma(1/x, shape1, scale = 1/scale).
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004 and
### <http://en.wikipedia.org/wiki/Inverse-gamma_distribution>
###
### AUTHORS:  Mathieu Pigeon, Christophe Dutang and
### Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvgamma <- function (x, shape, rate = 1, scale = 1/rate, log = FALSE)
     .External("actuar_do_dpq", "dinvgamma", x, shape, scale, log)

pinvgamma <- function(q, shape, rate = 1, scale = 1/rate,
                      lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pinvgamma", q, shape, scale, lower.tail, log.p)

qinvgamma <- function(p, shape, rate = 1, scale = 1/rate,
                      lower.tail = TRUE, log.p = FALSE)
     .External("actuar_do_dpq", "qinvgamma", p, shape, scale, lower.tail, log.p)

rinvgamma <- function(n, shape, rate = 1, scale = 1/rate)
     .External("actuar_do_random", "rinvgamma", n, shape, scale)

minvgamma <- function(order, shape, rate = 1, scale = 1/rate)
     .External("actuar_do_dpq", "minvgamma", order, shape, scale, FALSE)

levinvgamma <- function(limit, shape, rate = 1, scale = 1/rate,
                        order = 1)
     .External("actuar_do_dpq", "levinvgamma", limit, shape, scale, order, FALSE)

mgfinvgamma <- function(x, shape, rate = 1, scale = 1/rate, log = FALSE)
    .External("actuar_do_dpq", "mgfinvgamma", x, shape, scale, log)
