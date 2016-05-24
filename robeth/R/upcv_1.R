"Upcv" <-
function(svals)
{
        if(missing(svals)) return(6)
        n <- length(svals)
        fvals <- double(n)
        f.res <- .Fortran("upcva",
                n = to.integer(n),
                svals = to.single(svals),
                fvals = as.double(fvals))
        f.res$fvals
}

