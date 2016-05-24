### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {m,lev}beta functions to compute raw and limited
### moments for the Beta distribution (as defined in R). The
### noncentral beta distribution is _not_ supported.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

mbeta <- function(order, shape1, shape2)
    .External("actuar_do_dpq", "mbeta", order, shape1, shape2, FALSE)

levbeta <- function(limit, shape1, shape2, order = 1)
    .External("actuar_do_dpq", "levbeta", limit, shape1, shape2, order, FALSE)
