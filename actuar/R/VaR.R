### ===== actuar: An R Package for Actuarial Science =====
###
### Definition of a generic function for the Value at Risk. Currently,
### there is no default method for this function. The only method is
### for objects of class 'aggregateDist'.
###
### AUTHORS: Tommy Ouellet, Vincent Goulet <vincent.goulet@act.ulaval.ca>

VaR <- function(x, ...)
    UseMethod("VaR")

VaR.aggregateDist <- function(x, conf.level = c(0.9, 0.95, 0.99),
                              smooth = FALSE, names = TRUE, ...)
    quantile.aggregateDist(x, conf.level, smooth, names, ...)
