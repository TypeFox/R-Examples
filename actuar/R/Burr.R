### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}burr functions to compute
### characteristics of the Burr distribution. The version used in
### these functions has cumulative distribution function
###
###   Pr[X <= x] = 1 - (1/(1 + (x/scale)^shape2))^shape1, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dburr <- function (x, shape1, shape2, rate = 1, scale = 1/rate,
                   log = FALSE)
    .External("actuar_do_dpq", "dburr", x, shape1, shape2, scale, log)

pburr <- function(q, shape1, shape2, rate = 1, scale = 1/rate,
                  lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "pburr", q, shape1, shape2, scale,
              lower.tail, log.p)

qburr <- function(p, shape1, shape2, rate = 1, scale = 1/rate,
                  lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qburr", p, shape1, shape2, scale,
              lower.tail, log.p)

rburr <- function(n, shape1, shape2, rate = 1, scale = 1/rate)
    .External("actuar_do_random", "rburr", n, shape1, shape2, scale)

mburr <- function(order, shape1, shape2, rate = 1, scale = 1/rate)
    .External("actuar_do_dpq", "mburr", order, shape1, shape2, scale, FALSE)

levburr <- function(limit, shape1, shape2, rate = 1, scale = 1/rate,
                    order = 1)
    .External("actuar_do_dpq", "levburr", limit, shape1, shape2, scale, order, FALSE)
