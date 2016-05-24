"Psi" <-
function(svals)
{
        if(missing(svals)) return(1)
        n <- length(svals)
        fvals <- single(n)
        f.res <- .Fortran("psia",
                n = to.integer(n),
                svals = to.single(svals),
                fvals = to.single(fvals))
        f.res$fvals
}

