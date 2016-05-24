sospecden <- function (x, l, kernel = c("Trap", "Rect", "SupSm"), x.points = seq(-pi, pi, len = 200)) {
	kernel <- match.arg(kernel, c("Trap", "Rect", "SupSm"))
    kappa <- switch(kernel, Trap = kappaTrap, Rect = kappaRect, SupSm = kappaInDf)
    if (missing(l)) 
        l <- bwadap.ts(x)

	xacf <- as.vector(acf(x, lag.max = max((2 * l - 1), 0), type = "cov", 
        plot = F)$acf)

	Fones <- function(s) {
		M <- bw2order(xacf*kappa(0:(2*l-1), l), length(x))
		M <- min(M, length(x))
		xacf <- as.vector(acf(x, max(lag.max=(M-1), 0), type="cov", plot=F)$acf)
		(2*sum(kappaParzen(1:(M-1), M) * xacf[-1] * cos(s * 1:(M-1))) + xacf[1])/(2*pi)   
	}
	
	f <- Vectorize(Fones, "s")
    if (is.null(x.points)) {
        return(f)
    }
    return(list(x = x.points, y = f(x.points)))
}