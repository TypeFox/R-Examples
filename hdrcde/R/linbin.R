# Extracted from KernSmooth KernSmooth/R/all.R
## original file Copyright (C) M. P. Wand
## modifications for use with R copyright (C) B. D. Ripley
## Unlimited use and distribution (see LICENCE).

linbin <- function (X, gpoints, truncate = TRUE) 
{
    n <- length(X)
    M <- length(gpoints)
    if (truncate) 
        trun <- 1L
    else 0L
    a <- gpoints[1L]
    b <- gpoints[M]
    .Fortran("linbin", as.double(X), as.integer(n), as.double(a), 
        as.double(b), as.integer(M), as.integer(trun), double(M), PACKAGE="hdrcde")[[7]]
}
