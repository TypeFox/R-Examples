### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}invweibull functions to compute
### characteristics of the Inverse Weibull distribution. The version
### used in these functions has cumulative distribution function
###
###   Pr[X <= x] = exp(-(x/scale)^shape), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvweibull <- function (x, shape, rate = 1, scale = 1/rate, log = FALSE)
    .External("actuar_do_dpq", "dinvweibull", x, shape, scale, log)

pinvweibull <- function(q, shape, rate = 1, scale = 1/rate,
                        lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pinvweibull", q, shape, scale, lower.tail, log.p)

qinvweibull <- function(p, shape, rate = 1, scale = 1/rate,
                        lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qinvweibull", p, shape, scale, lower.tail, log.p)

rinvweibull <- function(n, shape, rate = 1, scale = 1/rate)
    .External("actuar_do_random", "rinvweibull", n, shape, scale)

minvweibull <- function(order, shape, rate = 1, scale = 1/rate)
    .External("actuar_do_dpq", "minvweibull", order, shape, scale, FALSE)

levinvweibull <- function(limit, shape, rate = 1, scale = 1/rate,
                          order = 1)
    .External("actuar_do_dpq", "levinvweibull", limit, shape, scale, order, FALSE)

## Aliases
dlgompertz <- dinvweibull
plgompertz <- pinvweibull
qlgompertz <- qinvweibull
rlgompertz <- rinvweibull
mlgompertz <- minvweibull
levlgompertz <- levinvweibull
