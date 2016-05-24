"Psp" <-
function(svals)
{
        if(missing(svals)) return(3)
        n <- length(svals)
        fvals <- single(n)
        f.res <- .Fortran("pspa",
                n = to.integer(n),
                svals = to.single(svals),
                fvals = to.single(fvals))
        f.res$fvals
}

