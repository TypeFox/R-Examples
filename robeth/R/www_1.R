"Www" <-
function(svals)
{
        if(missing(svals)) return(11)
        n <- length(svals)
        fvals <- double(n)
        f.res <- .Fortran("wwwa",
                n = to.integer(n),
                svals = to.single(svals),
                fvals = as.double(fvals))
        f.res$fvals
}

