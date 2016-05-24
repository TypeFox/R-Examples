"Chi" <-
function(svals)
{
        if(missing(svals)) return(4)
        n <- length(svals)
        fvals <- single(n)
        f.res <- .Fortran("chia",
                n = to.integer(n),
                svals = to.single(svals),
                fvals = to.single(fvals))
        f.res$fvals
}

