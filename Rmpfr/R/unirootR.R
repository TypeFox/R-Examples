### This is a translation of R_zeroin2 in ~/R/D/r-devel/R/src/appl/zeroin.c
### from  C to R  by John Nash,
### ---> file rootoned/R/zeroin.R of the new (2011-08-18) R-forge package rootoned
###
### Where John Nash calls it zeroin(), I call it unirootR()

##' Simple modification of uniroot()  which should work with mpfr-numbers
##' MM: uniroot() is in ~/R/D/r-devel/R/src/library/stats/R/nlm.R
##'
unirootR <- function(f, interval, ...,
		     lower = min(interval), upper = max(interval),
		     f.lower = f(lower, ...), f.upper = f(upper, ...),
		     verbose = FALSE,
		     tol = .Machine$double.eps^0.25, maxiter = 1000,
		     epsC = NULL)
{
    if(!missing(interval) && length(interval) != 2L)
	stop("'interval' must be a vector of length 2")
    ## For many "quick things", we will use as.numeric(.) but we do *NOT* assume that
    ## lower and upper are numeric!
    .N <- as.numeric
    if(.N(lower) >= .N(upper))
	##if(!is.numeric(lower) || !is.numeric(upper) || lower >= upper)
	stop("lower < upper  is not fulfilled")
    if(is.na(.N(f.lower))) stop("f.lower = f(lower) is NA")
    if(is.na(.N(f.upper))) stop("f.upper = f(upper) is NA")

    if((ff <- f.lower * f.upper) >= 0) {
	if(ff > 0)
	    stop("f() values at end points not of opposite sign")
	## else ff == 0	 <==> (at least) one of them is 0
	if(f.lower == 0)
	    return(list(root = lower, f.root = f.lower, iter = 0, estim.prec = tol))
	## else f.upper == 0 :
	return(list(root = upper, f.root = f.upper, iter = 0, estim.prec = tol))
    }

    if(is.null(epsC) || is.na(epsC)) {
	## determine 'epsC'  ``the achievable Machine precision''
	## -- given the class of f.lower, f.upper
	if(is.double(ff)) epsC <- .Machine$double.eps
	else if(is(ff, "mpfr"))
	    epsC <- 2^-min(getPrec(f.lower), getPrec(f.upper))
	else { ## another number class -- try to see if getPrec() is defined..
	    ## if not, there's not much we can do
	    if(is(prec <- tryCatch(min(getPrec(f.lower), getPrec(f.upper)),
				   error = function(e)e),
		  "error")) {
		warning("no valid getPrec() for the number class(es) ",
			paste(unique(class(f.lower),class(f.upper)), collapse=", "),
			".\n Using double precision .Machine$double.eps.")
		epsC <- .Machine$double.eps
	    } else {
		epsC <- 2^-prec
		message("using epsC = %s ..", format(epsC))
	    }
	}
    }
    if(tol < epsC)
        warning(sprintf("tol (%g) < epsC (%g)  is rarely sensical,
 and the resulting precision is probably not better than epsC",
                        tol, epsC))

    ## Instead of the call to C code, now "do it in R" :
    ## val <- .Internal(zeroin2(function(arg) as.numeric(f(arg, ...)),
    ##			     lower, upper, f.lower, f.upper,
    ##			     tol, as.integer(maxiter)))
    a <- lower # interval[1]
    b <- upper # interval[2]
    fa <- f.lower # f(ax, ...)
    fb <- f.upper # f(bx, ...)
    if (verbose)
	cat(sprintf("Start zeroin: f(%g)= %g;  f(%g)= %g\n",
		    .N(a), .N(fa),  .N(b), .N(fb)))
    c <- a
    fc <- fa
    ## First test if we have found a root at an endpoint
    maxit <- maxiter + 2 # count evaluations as maxiter-maxit
    converged <- FALSE

    while(!converged && maxit > 0) { ##---- Main iteration loop ------------------------------

	if (verbose) cat("Top of iteration, maxit=",maxit,"\n")
	d.prev <- b-a
	## Distance from the last but one to the last approximation	*/
	##double tol.2;		       ## Actual tolerance		*/
	##double p;			## Interpolation step is calcu- */
	##double q;			## lated in the form p/q; divi-
	##				 * sion operations is delayed
	##				 * until the last moment	*/
	##double d.new;		## Step at this iteration	*/

	if(abs(fc) < abs(fb)) { ## Swap data for b to be the smaller
	    if (verbose) cat("fc smaller than fb\n")
	    a <- b
	    b <- c
	    c <- a	## best approximation
	    fa <- fb
	    fb <- fc
	    fc <- fa
	}
	tol.2 <- 2*epsC*abs(b) + tol/2
	if (verbose) cat("tol.2= ",.N(tol.2),"\n")
	d.new <- (c-b)/2 # bisection
	if (verbose) cat("d.new= ",.N(d.new),"\n")

	converged <- abs(d.new) <= tol.2  ||  fb == 0
	if(converged) {
	    if (verbose) cat("DONE! -- small d.new or fb=0\n")
	    ## Acceptable approx. is found :
	    val <- list(root=b, froot=fb, rtol = abs(c-b), maxit=maxiter-maxit)
	}
	else {
	    ## Decide if the interpolation can be tried	*/
	    if( (abs(d.prev) >= tol.2) ## If d.prev was large enough*/
	       && (abs(fa) > abs(fb)) ) { ## and was in true direction,
		## Interpolation may be tried	*/
		##    register double t1,cb,t2;
		if (verbose) cat("d.prev larger than tol.2 and fa bigger than fb\n")

		cb <- c-b
		if (a == c) { ## If we have only two distinct points, linear interpolation
		    ## can only be applied
		    t1 <- fb/fa
		    p <- cb*t1
		    q <- 1 - t1
		    if (verbose) cat("a == c\n")
		}
		else { ## Quadric inverse interpolation*/
		    if (verbose) cat("a != c\n")
		    q <- fa/fc
		    t1 <- fb/fc
		    t2 <- fb/fa
		    p <- t2 * ( cb*q*(q-t1) - (b-a)*(t1-1) )
		    q <- (q-1) * (t1-1) * (t2-1)
		}
		if(p > 0) { ## p was calculated with the */
		    if (verbose) cat(" p > 0; ")
		    q <- -q ## opposite sign; make p positive */
		} else {    ## and assign possible minus to	*/
		    if (verbose) cat(" p <= 0; ")
		    p <- -p ## q				*/
		}
		if (p < 0.75*cb*q - abs(tol.2*q)/2 ## If b+p/q falls in [b,c]*/
		    && p < abs(d.prev*q/2)) {	## and isn't too large	*/
		    if (verbose) cat("p satisfies conditions for changing d.new\n")
		    d.new <- p/q ## it is accepted
		}
		## If p/q is too large, then the
		## bisection procedure can reduce [b,c] range to more extent
	    }

	    if( abs(d.new) < tol.2) { ## Adjust the step to be not less*/
		if (verbose) cat("d.new smaller than tol.2\n")
		if( d.new > 0 ) ## than tolerance		*/
		    d.new <- tol.2
		else
		    d.new <- -tol.2
	    }
	    a <- b
	    fa <- fb ## Save the previous approx. */
	    b <- b + d.new
	    fb <- f(b, ...)
	    maxit <- maxit-1
	    ## Do step to a new approxim. */
	    if( ((fb > 0) && (fc > 0)) || ((fb < 0) && (fc < 0)) ) {
		if (verbose) cat("make c to have sign opposite to b\n")
		## Adjust c for it to have a sign opposite to that of b */
		c <- a
		fc <- fa
	    }

	}## else

    } ## end{ while(maxit > 0) } --------------------------------------------

    if(converged) {
	iter <- val[["maxit"]]
	if(!is.na(fb) &&  abs(fb) > 0.5*max(abs(f.lower), abs(f.upper)))# from John Nash:
	    warning("Final function magnitude seems large -- maybe converged to sign-changing 'pole' location?")
    } else { ## (!converged) : failed!
	val <- list(root= b, rtol = abs(c-b))
	iter <- maxiter
	warning("_NOT_ converged in ", iter, " iterations")
    }

    list(root = val[["root"]], f.root = f(val[["root"]], ...),
	 iter = iter, estim.prec = .N(val[["rtol"]]), converged = converged)
}
