ed <- function (x, q = 1, w = 1, retq = TRUE)
{
    bsums <- function(x, q = 1, w = 1) {
        if (all(x == 0))
            return(0)
        if (q == 0)
            sum(x > 0)
        else if (q == 1)
            sum(-w * x * log(x), na.rm = TRUE)
        else sum(w * x^q)
    }
    rs <- rowSums(x)
	if (length(w) != 1 & length(w) != nrow(x))
		cat("Warning: n weights NE n rows of x !!")
    if (any(rs == 0)) {
        drops <- (1:nrow(x))[rs == 0]
        rs <- rs[rs != 0]
		x <- x[rs != 0,]
		w <- w[rs != 0]
	    cat("Dropping zero sum rows: ", drops, "\n")
    }
    x <- x/rs
    wa <- w^q/mean(w^q)
    a <- bsums(x, q = q, w = wa)/nrow(x)
    wg <- w/mean(w)
    g <- bsums(colMeans(wg * x, na.rm = TRUE), q = q, w = 1)
    if (retq) {
        if (q == 1) {
            a <- exp(a)
            g <- exp(g)
        }
        else if (q!= 1) {
            a <- a^(1/(1 - q))
            g <- g^(1/(1 - q))
        }
    c(alpha = a, beta = g/a, gamma = g)
    }
	else c(absums=a, gbsums=g)
}

