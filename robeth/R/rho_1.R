"Rho" <-
function(svals)
{
        if(missing(svals)) return(2)
        n <- length(svals)
        fvals <- single(n)
        f.res <- .Fortran("rhoa",
                n = to.integer(n),
                svals = to.single(svals),
                fvals = to.single(fvals))
        f.res$fvals
}

