"Ucv" <-
function(svals)
{
        if(missing(svals)) return(5)
        n <- length(svals)
        fvals <- double(n)
        f.res <- .Fortran("ucva",
                n = to.integer(n),
                svals = to.single(svals),
                fvals = as.double(fvals))
        f.res$fvals
}

