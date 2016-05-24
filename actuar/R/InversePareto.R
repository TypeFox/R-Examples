### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m}invpareto functions to compute
### characteristics of the Inverse Pareto distribution. The version
### used in these functions has cumulative distribution function
###
###   Pr[X <= x] = (x/(x + scale))^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dinvpareto <- function(x, shape, scale, log = FALSE)
     .External("actuar_do_dpq", "dinvpareto", x, shape, scale, log)

pinvpareto <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pinvpareto", q, shape, scale, lower.tail, log.p)

qinvpareto <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
     .External("actuar_do_dpq", "qinvpareto", p, shape, scale, lower.tail, log.p)

rinvpareto <- function(n, shape, scale)
     .External("actuar_do_random", "rinvpareto", n, shape, scale)

minvpareto <- function(order, shape, scale)
     .External("actuar_do_dpq", "minvpareto", order, shape, scale, FALSE)

levinvpareto <- function(limit, shape, scale, order = 1)
     .External("actuar_do_dpq", "levinvpareto", limit, shape, scale, order, FALSE)
