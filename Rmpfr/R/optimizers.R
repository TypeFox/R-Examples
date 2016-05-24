## From: Hans W Borchers <hwborchers@googlemail.com>
## To: Martin Maechler <maechler@stat.math.ethz.ch>
## Subject: optimizeR for Rmpfr
## Date: Sun, 3 Jun 2012 16:58:12 +0200

## This is from Hans' pracma package,
## /usr/local/app/R/R_local/src/pracma/R/golden_ratio.R
## but there's also fibonacci search,  direct1d, ....
optimizeR <- function(f, lower, upper, ..., tol = 1e-20,
		      method = c("Brent", "GoldenRatio"),
		      maximum = FALSE,
		      precFactor = 2.0,
		      precBits = -log2(tol) * precFactor, maxiter = 1000,
		      trace = FALSE)
{
    stopifnot(length(lower) == 1, length(upper) == 1, lower <= upper)
    fun <- match.fun(f)
    f <- if(maximum) function(x) -fun(x, ...) else function(x) fun(x, ...)

    a <- if(!is(lower,"mpfr")) mpfr(lower, precBits = precBits)
	 else if(.getPrec(lower) < precBits) roundMpfr(lower, precBits)
    b <- if(!is(upper,"mpfr")) mpfr(upper, precBits = precBits)
	 else if(.getPrec(upper) < precBits) roundMpfr(upper, precBits)
    method <- match.arg(method)
    n <- 0; convergence <- TRUE
    ## if(method == "GoldenRatio") {
    switch(method,
	   "GoldenRatio" =
       { ## golden ratio
	   phi <- 1 - (sqrt(mpfr(5, precBits = precBits)) - 1)/2
	   x <- c(a, a + phi*(b-a), b - phi*(b-a), b)
	   y2 <- f(x[2])
	   y3 <- f(x[3])
	   while ((d.x <- x[3] - x[2]) > tol) {
	       n <- n + 1
	       if(trace && n %% trace == 0)
		   message(sprintf("it.:%4d, delta(x) = %12.8g", n, as.numeric(d.x)))
	       if (y3 > y2) {
		   x[2:4] <- c(x[1]+phi*(x[3]-x[1]), x[2:3])
		   y3 <- y2
		   y2 <- f(x[2])
	       } else {
		   x[1:3] <- c(x[2:3], x[4]-phi*(x[4]-x[2]))
		   y2 <- y3
		   y3 <- f(x[3])
	       }
	       if (n > maxiter) {
		   warning(sprintf("not converged in %d iterations (d.x = %g)",
				   maxiter, as.numeric(d.x)))
		   convergence <- FALSE
		   break
	       }
	   }
	   xm <- (x[2]+x[3])/2
	   fxm <- if (abs(f. <- f(xm)) <= tol^2) 0. else f.
       },
	   "Brent" =
       {
	   ##--- Pure R version (for "Rmpfr") of   R's fmin() C code.

    ##	       The method used is a combination of  golden  section  search  and
    ##	   successive parabolic interpolation.	convergence is never much slower
    ##	   than	 that  for  a  Fibonacci search.  If  f	 has a continuous second
    ##	   derivative which is positive at the minimum (which is not  at  ax  or
    ##	   bx),	 then  convergence  is	superlinear, and usually of the order of
    ##	   about  1.324....

    ##	       The function  f	is never evaluated at two points closer together
    ##	   than	 eps*abs(fmin)+(tol/3), where eps is the square
    ##	   root	 of  2^-precBits.	  if   f   is a unimodal
    ##	   function and the computed values of	 f   are  always  unimodal  when
    ##	   separated  by  at least  eps*abs(x)+(tol/3), then  fmin  approximates
    ##	   the abcissa of the global minimum of	 f  on the interval  ax,bx  with
    ##	   an error less than  3*eps*abs(fmin)+tol.  if	  f   is  not  unimodal,
    ##	   then fmin may approximate a local, but perhaps non-global, minimum to
    ##	   the same accuracy.

    ##	       This function subprogram is a slightly modified	version	 of  the
    ##	   Algol  60 procedure	localmin  given in Richard Brent, Algorithms for
    ##	   Minimization without Derivatives, Prentice-Hall, Inc. (1973).

	   ## c is the squared inverse of the golden ratio
	   c <- (3 - sqrt(mpfr(5, precBits = precBits))) / 2

	   eps <- 2^-precBits
	   tol1 <- 1+eps		# the smallest 1.000... > 1
	   eps <- sqrt(eps)

	   w  <-  v <-	x <- a + c * (b - a)
	   fw <- fv <- fx <- f(x)
	   d <- e <- 0
	   tol3 <- tol / 3

	   ##  main loop starts here -----------------------------------
	   repeat {
	       n <- n+1
	       xm <- (a + b) /2
	       tol1 <- eps * abs(x) + tol3
	       t2 <- tol1 * 2

	       ## check stopping criterion
	       if (abs(x - xm) <= t2 - (d.x <- (b - a)/2))
		   break
	       if (n > maxiter) {
		   warning(sprintf("not converged in %d iterations (d.x = %g)",
				   maxiter, as.numeric(d.x)))
		   convergence <- FALSE
		   break
	       }
	       p <- q <- r <- 0
	       if (abs(e) > tol1) { ## fit parabola
		   r <- (x - w) * (fx - fv)
		   q <- (x - v) * (fx - fw)
		   p <- (x - v) * q - (x - w) * r
		   q <- (q - r) * 2
		   if (q > 0) p <- -p else q <- -q
		   r <- e; e <- d
	       }

	       if(doTr <- (trace && n %% trace == 0))
		   msg <- sprintf("it.:%4d, x = %-19.12g, delta(x) = %9.5g",
				  n, as.numeric(x), as.numeric(d.x))
	       if (abs(p) >= abs(q/2 * r) || p <= q * (a - x) ||
		   p >= q * (b - x)) { ## a golden-section step
		   e <- (if(x < xm) b else a) - x
		   d <- c * e
		   if(doTr) msg <- paste(msg, "+ Golden-Sect.")
	       }
	       else { ## a parabolic-interpolation step
		   d <- p / q
		   u <- x + d
		   if(doTr) msg <- paste(msg, "+ Parabolic")

		   ## f must not be evaluated too close to ax or bx
		   if (u - a < t2 || b - u < t2) {
		       d <- tol1
		       if (x >= xm) d <- -d
		   }
	       }
	       if(doTr) message(msg)

	       ## f must not be evaluated too close to x
	       u <- x + if(abs(d) >= tol1) d else if(d > 0) tol1 else -tol1
	       fu <- f(u)

	       ##  update  a, b, v, w, and x
	       if (fu <= fx) {
		   if (u < x) b <- x else a <- x
		   v <- w;    w <- x;	x <- u
		   fv <- fw; fw <- fx; fx <- fu
	       } else {
		   if (u < x) a <- u else b <- u
		   if (fu <= fw || w == x) {
		       v <- w; fv <- fw
		       w <- u; fw <- fu
		   } else if (fu <= fv || v == x || v == w) {
		       v <- u; fv <- fu
		   }
	       }
	   } ## end {repeat} main loop

	   xm <- x; fxm <- fx
       }, ## end{ "Brent" }

	   stop(sprintf("Method '%s' is not implemented (yet)", method)))

    c(if(maximum)list(maximum = xm) else list(minimum = xm),
      list(objective = fxm, iter = n,
	   convergence = convergence, estim.prec = abs(d.x), method=method))
}
