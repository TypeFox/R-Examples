##
##  r i d d e r s . R  Ridders' Method
##


ridders <- function(f, a, b, maxiter = 100, tol = .Machine$double.eps^0.5) {
    # (!is.numeric(a) && !is.complex(a) && !is(a,"mpfr") ||
    #  !is.numeric(b) && !is.complex(b) && !is(b,"mpfr"))
    #     stop("Arguments 'a' and 'b' must be numeric, complex, or mpfr.")
	x1 <- a;      x2 <- b
	f1 <- f(x1);  f2 <- f(x2)
	if (f1*f2 >= 0) stop("f(a) and f(b) must have different signs.")

    niter <- 2
    while(abs(x1 - x2) > tol && niter < maxiter) {
        xm <- (x1 + x2)/2;  fm <- f(xm)
        if (fm == 0)
            return(list(root = xm, f.root = 0, niter = niter, estim.prec = 0))

        x3 <- xm + (xm - x1) * sign(f1 - f2) * fm / sqrt(fm^2 - f1 * f2)
        f3 <- f(x3);  niter <- niter + 2
        if (f3 == 0)
            return(list(root = x3, f.root = 0, niter = niter, estim.prec = 0))

        if (fm * f3 < 0) {
            x1 <- xm;  f1 <- fm
            x2 <- x3;  f2 <- f3
        } else if (f1 * f3 < 0) {
            x2 <- x3;  f2 <- f3
        } else if (f2 * f3 < 0) {
            x1 <- x3;  f1 <- f3
        } else {
            stop("Inform the maintainer: you should never get here.")
        }
     }

    if (abs(f1) < abs(f2)) {
        x0 <- x1;  f0 <- f1
    } else {
        x0 <- x2;  f0 <- f2
    }
    ep <- abs(x1 - x2)
    return(list(root = x0, f.root = f0, niter = niter, estim.prec = ep))
}
