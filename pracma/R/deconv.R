##
##  d e c o n v . R  Deconvolution
##


deconv <- function(b, a) {
    if (length(b) == 0)
        return(list(q = 0, r = c()))
    if ( (!is.numeric(b) && ! is.complex(b)) ||
         (!is.numeric(a) && ! is.complex(a)) )
        stop("Arguments 'b' and 'a' must be numeric or complex.")

    if ( a[1] == 0)
        stop("First element of argument 'a' must be nonzero.")

    nb <- length(b)
    na <- length(a)
    if (nb < na)
        return(list(q = 0, r = b))

    q <- c()
    while (nb >= na) {
        d <- b[1] / a[1]
        b <- b - conv(a, c(d, rep(0, nb-na)))
        q <- c(q, d)
        b <- b[2:nb]
        nb <- nb -1
    }
    return(list(q = q, r = b))
}
