maximType <- function(x)
    UseMethod("maximType")

maximType.default <- function(x)
    x$maximType

maximType.maxim <- function(x)
    x$type
