### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}invexp functions to compute
### characteristics of the Inverse Exponential distribution. The
### version used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = exp(-scale/x), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvexp <- function (x, rate = 1, scale = 1/rate, log = FALSE)
    .External("actuar_do_dpq", "dinvexp", x, scale, log)

pinvexp <- function(q, rate = 1, scale = 1/rate,
                    lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pinvexp", q, scale, lower.tail, log.p)

qinvexp <- function(p, rate = 1, scale = 1/rate,
                    lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qinvexp", p, scale, lower.tail, log.p)

rinvexp <- function(n, rate = 1, scale = 1/rate)
    .External("actuar_do_random", "rinvexp", n, scale)

minvexp <- function(order, rate = 1, scale = 1/rate)
    .External("actuar_do_dpq", "minvexp", order, scale, FALSE)

levinvexp <- function(limit, rate = 1, scale = 1/rate, order)
    .External("actuar_do_dpq", "levinvexp", limit, scale, order, FALSE)
