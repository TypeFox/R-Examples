# Simply call 'first' or 'last' with a different default value for 'n'.
first <- function(x, n=1, ...) head(x, n=n, ...)
last  <- function(x, n=1, ...) tail(x, n=n, ...)
