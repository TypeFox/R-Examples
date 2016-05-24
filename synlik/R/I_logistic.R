# Logistic                             R -> (0, 1)
# and inverse logistic (logit)         (0, 1) -> R 
# transformations

.logistic  <- function(x) 1 / (1 + exp(-x))

.ilogistic <- function(x) log( x / (1 - x) )