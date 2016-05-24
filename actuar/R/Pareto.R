### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {d,p,q,r,m,lev}pareto functions to compute
### characteristics of the Pareto distribution. The version used in
### these functions has cumulative distribution function
###
###   Pr[X <= x] = 1 - (scale/(x + scale))^shape, x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

dpareto <- function (x, shape, scale, log = FALSE)
    .External("actuar_do_dpq", "dpareto", x, shape, scale, log)

ppareto <- function (q, shape, scale, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "ppareto", q, shape, scale, lower.tail, log.p)

qpareto <- function (p, shape, scale, lower.tail = TRUE, log.p = FALSE)
    .External("actuar_do_dpq", "qpareto", p, shape, scale, lower.tail, log.p)

rpareto <- function(n, shape, scale)
    .External("actuar_do_random", "rpareto", n, shape, scale)

mpareto <- function(order, shape, scale)
     .External("actuar_do_dpq", "mpareto", order, shape, scale, FALSE)

levpareto <- function(limit, shape, scale, order = 1)
     .External("actuar_do_dpq", "levpareto", limit, shape, scale, order, FALSE)

## Aliases
dpareto2 <- dpareto
ppareto2 <- ppareto
qpareto2 <- qpareto
rpareto2 <- rpareto
mpareto2 <- mpareto
levpareto2 <- levpareto
