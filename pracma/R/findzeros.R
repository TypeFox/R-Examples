##
##  f i n d z e r o s . R  Find all roots or minima
##


findzeros <- function(f, a, b, n = 100, tol = .Machine$double.eps^(2/3), ...) {
    stopifnot(is.numeric(a), length(a) == 1,
              is.numeric(b), length(b) == 1,
              is.numeric(n), floor(n) == ceiling(n), n >= 2)
    if (! a < b)
        stop("Left interval border must be smaller than right one.")

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

	h <- (b - a) / n
	x <- seq(a, b, by = h)  # length(x) == n+1
	y <- f(x)

    R <- c()
    s <- sign(f(x[1]))
    if (abs(f(x[1])) < tol) {
        R <- c(x[1])
        s <- 0
    }

	for (i in 2:n) {
	    si <- sign(f(x[i]))
	    if (abs(f(x[i])) < tol) {
	        R <- c(R, x[i])
	        si <- 0
	    } else if (s * si < 0) {  # function values have different sign, != 0
		    u <- uniroot(f, c(x[i-1], x[i]))
		    R <- c(R, u$root)
		} else if (s * si > 0) {  # function values both positive or negative
		    xm <- (x[i-1] + x[i])/2
		    ym <- f(xm)
		    d <- (y[i] - y[i-1])/h
		    if (d == 0) next
		    xv <- xm - ym/d
		    if (xv > x[i-1] && xv < x[i]) {
		        if (s > 0) {
		            s <- optimize(f, c(x[i-1], x[i]), tol = tol)
		            sm <- s$minimum
	            } else {
	                s <- optimize(f, c(x[i-1], x[i]), maximum = TRUE, tol = tol)
	                sm <- s$maximum
	            }
		        if (abs(s$objective) < tol)
		            R <- c(R, sm)
		    }
		}
		s <- si
	}
    if (abs(f(x[n+1])) < tol)
        R <- c(R, x[n+1])

	return(R)
}


findmins <- function(f, a, b, n = 100, tol = .Machine$double.eps^(2/3), ...) {
    stopifnot(is.numeric(a), length(a) == 1,
              is.numeric(b), length(b) == 1,
              is.numeric(n), floor(n) == ceiling(n), n >= 2)
    if (! a < b)
        stop("Left interval border must be smaller than right one.")

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

	h <- (b - a) / n
	x <- seq(a, b, by = h)  # length(x) == n+1

    R <- c()

    for (i in 2:(n-1)) {
        if ( (f(x[i]) - f(x[i-1]) < 0) && (f(x[i+1]) - f(x[i])) > 0 ) {
            o <- optimize(f, c(x[i-1], x[i+1]))
            R <- c(R, o$minimum)
        }
    }

    return(R)
}
