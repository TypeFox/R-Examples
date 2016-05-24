eds <- function(x, q = 1, w = 1, retq = TRUE) {
	bsums <- function(x, w) {
		if (all(x == 0)) return(0)
		sum(-w*x*log(x), na.rm=TRUE)
	}
	if (is.null(dim(x))) x <- matrix(x,nrow=1)
	x[x!=0] <- x[x!=0]^q
	rs <- rowSums(x)
	if (length(w) != 1 & length(w) != nrow(x))
		cat("Warning: n weights NE n rows of x !!")
	if (any(rs==0)) {
		drops <- (1:nrow(x))[rs==0]
		rs <- rs[rs!=0]
		w <- w[rs != 0]
		x <- x[rs!=0,]
		cat("Dropping zero sum rows: ",drops,"\n")
	}
	x <- x/rs
    wa <- w^q/mean(w^q)
    a <- bsums(x, w = wa)/nrow(x)
    wg <- w/mean(w)
    g <- bsums(colMeans(wg * x, na.rm = TRUE), w = 1)
	if (retq) {
		a <- exp(a)
		g <- exp(g)
		c(alpha=a, beta=g/a, gamma=g)
	}
	else c(absums=a, gbsums=g)
}

