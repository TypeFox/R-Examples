"Wcv" <-
function(svals)
{
        if(missing(svals)) return(9)
        n <- length(svals)
        fvals <- double(n)
        f.res <- .Fortran("wcva",
                n = to.integer(n),
                svals = to.single(svals),
                fvals = as.double(fvals))
        f.res$fvals
}

