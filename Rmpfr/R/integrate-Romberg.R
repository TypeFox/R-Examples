####  Romberg integration in pure R
####  ===================    ====== so it can be used with  Rmpfr

## TODO: Lauren K  would like to get return()ed all the intermediate sums as well
## ----  I agree: but only if  'all.sums = TRUE'

integrateR <- function(f, lower, upper, ..., ord = NULL,
		       rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
		       max.ord = 19, verbose = FALSE)
{
    stopifnot(length(lower) == 1, length(upper) == 1,
              is.finite(lower), is.finite(upper))
    f <- match.fun(f)
    ff <-
	## if(verbose) function(x) { cat("f(x), x="); str(x) ; f(x, ...) } else
	function(x) f(x, ...)
    null.ord <- is.null(ord) ## ord := Romberg order
    has.tol <- !missing(rel.tol) || !missing(abs.tol)# one of them specified
    if(chkConv <- (null.ord || has.tol)) {
	## will use rel.tol and abs.tol
	if (abs.tol <= 0 && rel.tol < max(50 * .Machine$double.eps, 5e-29))
	    stop("invalid tolerance values")
	## but need (maximal) order for Bauer's algorithm:
    }
    if(null.ord && !has.tol) {
	## This is "approximate" (and too large, typically); but if it's
	## too small, t[.] will be extended automatically:
	ord <- ## == max(3, min(25, ceiling(-log2(rel.tol))))  with default rel.tol
	    13
	if(verbose) cat("ord =", ord, "(as no '*.tol' specified)\n")
    }
    useOrd <- !null.ord || !has.tol
    if(useOrd) {
	stopifnot(ord >= 0)
	if(verbose) cat(sprintf(
	    " ord = %d; ==> evaluating integrand at %s 2^(ord+1)-2 = %d locations\n",
	    ord, if(chkConv) "up to" else "", 2^(ord+1)-2))
    }

### Bauer(1961)  "Algorithm 60 -- Romberg Integration" Comm.ACM 4(6), p.255
    m <- le <- upper - lower # 'l'
    ##  a "hack", but really improves the result:
    if(!is.numeric(m)) {
	if(is(m, "mpfr")) { # should get *same* precision
	    if(is.numeric(lower)) lower <- 0*m+ lower
	    if(is.numeric(upper)) upper <- 0*m+ upper
	} else { ## other high-precision...
	    if(is.numeric(lower)) lower <- as(lower, class(m))
	    if(is.numeric(upper)) upper <- as(upper, class(m))
	}
    }
    t1 <- (ff(lower) + ff(upper))/2
    t <- rep(t1, if(useOrd) ord+1 else 10)## <- must work for "mpfr" numbers
    one <- 1 + 0*t1 # for "mpfr"

    r. <- t[1]*le
    n <- 1 # 'n'(Bauer) = 2^n (Romberg Algo)
    if(verbose) { ## In "mpfr" case, cannot use sprintf("%g");
	## ==> rather use format(.) with higher number of digits
	prDigs <- max(10, min(50, 2 + ceiling(-log10(rel.tol))))
	FORM <- paste0("n=%2d, 2^n=%9.0f | I = %",(5+prDigs), "s, abs.err =%14s\n")
    }
    h <- 1L
    repeat {
	if(verbose >= 2) {
	    cat("range(le*t[1:h]):\n\t");
	    print(format(range(le*t[1:h]), digits=15), quote=FALSE)
	}
	u <- 0
	m <- m/2 # == le/(2*n)
	## here, we require f(.) to be vectorized:
	u <- sum(ff(lower+ seq(1, 2*n-1, by=2)*m))
	t[h+1L] <- (u/n + t[h])/2
	f. <- one
	for(j in h:1) {
	    f. <- 4*f.
	    t[j] <- t[j+1] + (t[j+1] - t[j])/ (f. - 1)
	}
        r <- t[1]*le
        aErr <- abs(r - r.)
	if(verbose)
	    cat(sprintf(FORM, h, 2*n, format(r, digits = prDigs),
			format(aErr, digits = max(7, getOption("digits")))))
	if(chkConv) { ## check if we converged: |r - r.| < min(.,.): *both* tolerances must be satisfied
	    if(converged <- (aErr < min(abs(r)*rel.tol, abs.tol)))
		break
	}
        if((useOrd && h >= ord) || h >= max.ord) break
        r. <- r
	n <- 2*n # == 2^h
        h <- h+1L
    }
    if(chkConv && !converged) {
	relE <- format(aErr/abs(r), digits=getOption("digits"))
	msg <- paste0("no convergence up to order ", ord,
		      "; last relative change = ", relE, "\n",
                      "Consider setting 'ord = <n>' (e.g. = ", ord+1,").")
	warning(msg)
    } else msg <- "OK"
    r <- list(value = r, abs.error = aErr, subdivisions = 2*n+1,
              "message" = msg, call = match.call())
    class(r) <- c("integrateR", "integrate")
    r
}

## This is  such that  print.integrate() calls our format(<mpfr>) method
## (and do *not* hide it via S3method() in NAMESPACE):
## print.integrate <- getS3method("print","integrate")# from 'stats' possibly not exported
## environment(print.integrate) <- environment()
## setMethod(show, "integrate", function(object) print.integrate(object))

print.integrateR <- function (x, digits = max(3, getOption("digits")-2), ...)
{
    if(x[["message"]] != "OK")
        cat("Non-convergence message ", sQuote(x$message), "\n", sep = "")
    ## The "Ok" message:
    cat(format(x$value, digits = digits),
       " with absolute error < ", format(x$abs.error, digits=digits),
       "\n", sep = "")
    invisible(x)
}
setMethod(show, "integrateR", function(object) print.integrateR(object))
