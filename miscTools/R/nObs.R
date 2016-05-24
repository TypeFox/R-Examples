## Return #of observations for models
nObs <- function(x, ...)
    ## Number of observations for statistical models
    UseMethod("nObs")

nObs.lm <- function(x, ...)
    nrow(x$qr$qr)

nObs.default <- function(x, ...)
    x$param$nObs

