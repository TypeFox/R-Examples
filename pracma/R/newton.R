##
##  n e w t o n . R  Newton Root finding
##


newtonRaphson <- function(fun, x0, dfun = NULL, ...,
	               maxiter = 100, tol = .Machine$double.eps^0.5) {
	# Newton method for finding function zeros
	if  (is.null(dfun)) {
		dfun <- function(x, ...) { h <- tol^(2/3)
			(fun(x+h, ...) - fun(x-h, ...)) / (2*h)
		}
	}
	x   <- x0
	fx  <- fun(x, ...)
	dfx <- dfun(x, ...)
	niter <- 0
	diff  <- tol + 1
	while (diff >= tol && niter <= maxiter) {
		niter <- niter + 1
        if (dfx == 0) {
            warning("Slope is zero: no further improvement possible.")
            break
        }
		diff  <- - fx/dfx
		x <- x + diff
		diff <- abs(diff)
		fx  <- fun(x, ...)
		dfx <- dfun(x, ...)
	}
	if (niter > maxiter) {
		warning("Maximum number of iterations 'maxiter' was reached.")
	}
	return(list(root=x, f.root=fx, niter=niter, estim.prec=diff))
}


newton <- newtonRaphson


halley <- function(fun, x0,
                   maxiter = 100, tol = .Machine$double.eps^0.5) {
    f0 <- fun(x0)
    if (abs(f0) < tol^(3/2))
        return(list(root = x0, f.root = f0, maxiter = 0, estim.prec = 0))

    f1 <- fderiv(fun, x0, 1)
    f2 <- fderiv(fun, x0, 2)
    x1 <- x0 - 2*f0*f1 / (2*f1^2 - f0*f2)

    niter = 1
    while (abs(x1 - x0) > tol && niter < maxiter) {
        x0 <- x1
        f0 <- fun(x0)
        f1 <- fderiv(fun, x0, 1)
        f2 <- fderiv(fun, x0, 2)
        x1 <- x0 - 2*f0*f1 / (2*f1^2 - f0*f2)
        niter <- niter + 1
    }
    return(list(root = x1, f.root = fun(x1),
                iter = niter, estim.prec = abs(x1 - x0)))
}


newtonHorner <- function(p, x0, 
                         maxiter = 50, tol = .Machine$double.eps^0.5) {
    n <- length(p) - 1
    niter <- 0
    x <- x0
    diff <- 1 + tol

    while (niter <= maxiter && diff >= tol) {
        H <- horner(p, x)
        if (abs(H$dy) <= tol) {
            warning("Newton's method encountered a slope almost zero.")
            return(list(root = NULL, f.root = NULL, deflate = NULL,
                        iters = niter, estim.prec = Inf))
        }
        xnew <- x - H$y / H$dy
        diff <- abs(x - xnew)
        niter <- niter + 1
        x <- xnew
    }

    if (niter > maxiter) {
        warning("Maximum number of iterations exceeded.")
    }
    defl <- hornerdefl(p, x)

    return(list(root = x, f.root = defl$y, deflate = defl$q,
                iters = niter, estim.prec = diff))
}

