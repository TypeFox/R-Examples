gauss.quad <- ## Generate Gauss-Legendre quadrature
function(size,interval) {
    if (interval[1]>=interval[2])
        warning("gss warning in gauss.quad: interval limits swapped")
    z <- .Fortran("gaussq",
                  as.integer(1),
                  as.integer(size),
                  as.double(0), as.double(0),
                  as.integer(0),
                  as.double(c(-1,1)), double(size),
                  t=double(size), w=double(size),
                  PACKAGE="gss")
    mn <- min(interval[1:2])
    range <- abs(interval[1]-interval[2])
    pt <- mn+range*(z$t+1)/2
    wt <- range*z$w/2
    list(pt=pt,wt=wt)
}
