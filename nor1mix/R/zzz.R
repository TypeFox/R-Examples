## no S4 methodology here; speedup :
.noGenerics <- TRUE

if(FALSE) # no longer {was '.First.lib', .. and hence unused from ~ 2011}
.onLoad <- function(libname, pkg) {
    if(interactive() || getOption("verbose")) { # not in test scripts
	packageStartupMessage(sprintf("Package %s (%s) loaded.  To cite, see citation(\"%s\")\n",
			pkg, packageDescription(pkg)$Version, pkg))
## was end of 2007:
	warning("'MW.nm2' has been changed (linear transform only) to match\n",
		"the Annals paper. 'MW.nm2.old' is the former version",
		call. = FALSE, immediate. = TRUE)
    }
}


## For reference; not needed anymore, as ../DESCRIPTION now has  R >= 2.8.0 :
if(getRversion() < "2.8.0") {
monoHsplinefun <- function(x, y=NULL, ties = mean)
{
    ## Purpose: "Monotone cubic Hermite interpolation"
    ## 		like stats::splinefun()  but using *monotone* splines
    ## ----------------------------------------------------------------------
    ## Arguments: x,y: numeric (data)
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  9 Jan 2008, 15:10

    ##  *very much* like splinefun() ..../R/src/library/stats/R/splinefun.R :
    x <- xy.coords(x, y)
    y <- x$y
    x <- x$x
    n <- length(x)
    if(any(o <- is.na(x) | is.na(y))) {
	o <- !o
	x <- x[o]
	y <- y[o]
	n <- length(x)
    }
    if (!identical(ties, "ordered")) {
	if (length(ux <- unique(x)) < n) {
	    if (missing(ties))
		warning("collapsing to unique 'x' values")
	    y <- as.vector(tapply(y,x,ties))# as.v: drop dim & dimn.
	    x <- sort(ux)
	    n <- length(x)
	    rm(ux)
	} else {
	    o <- order(x)
	    x <- x[o]
	    y <- y[o]
	}
    }
    if(n == 0) stop("zero non-NA points")

    n1 <- n - 1L
    ## - - - "Data preprocessing" - - -

    dy <- y[-1] - y[-n] # = diff(y)
    i0 <- dy == 0 # or |dy| < eps ?? fixme ??
    dx <- x[-1] - x[-n] # = diff(x)
    Sx <- dy / dx       # 2. \Delta_k = (y_{k+1} - y_k)/(x_{k+1} - x_k), k=1:n1
    m <- c(Sx[1], (Sx[-1] + Sx[-n1])/2, Sx[n1]) ## 1.
    if(any(i0)) {
        ## m0[k] := i0[k] or i0[k-1]
        m0 <- c(i0,FALSE) | c(FALSE,i0)
        m[m0] <- 0
    }
    if(any(ip <- !i0)) {
        alpha <- m[-n][ip] / Sx[ip]
        beta  <- m[-1][ip] / Sx[ip]
        a2b3 <- 2*alpha + beta - 3
        ab23 <- alpha + 2*beta - 3
        if(any(ok <- (a2b3 > 0 & ab23 > 0)) &&
           any(ok <- ok & (alpha * (a2b3 + ab23) < a2b3^2))) {
            tau <- 3 / sqrt(alpha[ok]^2 + beta[ok]^2)
            m[-n][ip][ok] <- tau * alpha[ok] * Sx[ip][ok]
            m[-1][ip][ok] <- tau *  beta[ok] * Sx[ip][ok]
        }
    }
    ## now do "Hermite spline with (x,y,m)":
    Hsplinefun0(x = x, y = y, m = m, dx = dx)
}

## internal function; the user function is  Hsplinefun()
Hsplinefun0 <- function(x, y, m, dx = x[-1L] - x[-length(x)])
{
    ## return a
    function(u, deriv=0, extrapol = c("linear","cubic"))
    {
        extrapol <- match.arg(extrapol)
	deriv <- as.integer(deriv)
	if (deriv < 0 || deriv > 2)
	    stop("'deriv' must be between 0 and 2")
	i <- findInterval(u, x, all.inside = (extrapol == "cubic"))
	if(deriv == 0)
	    interp <- function(u, i) {
		h <- dx[i]
		t <- (u - x[i]) / h
		## Compute the 4 Hermite (cubic) polynomials h00, h01,h10, h11
		t1 <- t-1
		h01 <- t*t*(3 - 2*t)
		h00 <- 1 - h01
		tt1 <- t*t1
		h10 <- tt1 * t1
		h11 <- tt1 * t
		y[i]  * h00 + h*m[i]  * h10 +
		y[i+1]* h01 + h*m[i+1]* h11
	    }
	else if(deriv == 1)
	    interp <- function(u, i) {
		h <- dx[i]
		t <- (u - x[i]) / h
		## 1st derivative of Hermite polynomials h00, h01,h10, h11
		t1 <- t-1
		h01 <- -6*t*t1 # h00 = - h01
		h10 <- (3*t - 1) * t1
		h11 <- (3*t - 2) * t
		(y[i+1] - y[i])/h * h01 + m[i] * h10 + m[i+1]* h11
	    }
	else ## deriv == 2
	    interp <- function(u, i) {
		h <- dx[i]
		t <- (u - x[i]) / h
		## 2nd derivative of Hermite polynomials h00, h01,h10, h11
		h01 <- 6*(1-2*t) # h00 = - h01
		h10 <- 2*(3*t - 2)
		h11 <- 2*(3*t - 1)
		((y[i+1] - y[i])/h * h01 + m[i] * h10 + m[i+1]* h11) / h
	    }


	if(extrapol == "linear" &&
	   any(iXtra <- (iL <- (i == 0)) | (iR <- (i == (n <- length(x)))))) {
	    ##	do linear extrapolation
	    r <- u
	    if(any(iL)) r[iL] <- if(deriv == 0) y[1] + m[1]*(u[iL] - x[1]) else
				  if(deriv == 1) m[1] else 0
	    if(any(iR)) r[iR] <- if(deriv == 0) y[n] + m[n]*(u[iR] - x[n]) else
				  if(deriv == 1) m[n] else 0
	    ## For internal values, compute "as normal":
	    ini <- !iXtra
	    r[ini] <- interp(u[ini], i[ini])
	    r
	}
        else { ## use cubic Hermite polynomials, even for extrapolation
            interp(u, i)
        }

    }
}

}## R version < 2.8.0
