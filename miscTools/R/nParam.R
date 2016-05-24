## Return the #of parameters of model
nParam <- function(x, free=FALSE, ...)
    UseMethod("nParam")

nParam.default <- function(x, ...)
    x$param$nParam

nParam.lm <- function(x, ...)
    length(coefficients(x))

