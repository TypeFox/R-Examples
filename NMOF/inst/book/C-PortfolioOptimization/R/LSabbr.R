# LSabbr.R -- version 2010-12-31
function(OF, algo=list(), ...) {
    # add defaults to algo, do some checks
    # ...
    OF1 <- function(x) OF(x,...)
    N1 <- function(x) algo$neighbour(x,...)
    Fmat <- array(NA, dim = c(algo$nS,2))
    xc  <- algo$x0; xcF <- OF1(xc)
    for (s in 1:algo$nS){	
        xn <- N1(xc)
        xnF <- OF1(xn)
        if (xnF <= xcF){
            xc  <- xn
            xcF <- xnF
        }
        Fmat[s,1] <- xnF
        Fmat[s,2] <- xcF
    }
    return( list(xbest = xc, OFvalue = xcF, Fmat = Fmat) )
}