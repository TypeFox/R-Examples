### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of the {m,lev}weibull functions to compute raw and
### limited moments for the Weibull distribution (as defined in R).
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mweibull <- function(order, shape, scale = 1)
    .External("actuar_do_dpq", "mweibull", order, shape, scale, FALSE)

levweibull <- function(limit, shape, scale = 1, order = 1)
    .External("actuar_do_dpq", "levweibull", limit, shape, scale, order, FALSE)
