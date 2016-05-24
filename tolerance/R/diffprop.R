ddiffprop <- function(x, k1, k2, n1, n2, a1 = 0.5, a2 = 0.5, log = FALSE, ...)
{
    if ((a1 < 0 | a1 > 1) | (a2 < 0 | a2 > 1)) {
        warning("a1 and a2 should both be between 0 and 1 for this fiducial approach!", call. = TRUE)
    }
	d <- x
	ddiffprop2 <- function(d, k1, k2, n1, n2, a1, a2,...)
	{
		c1 <- k1 + a1
		c2 <- k2 + a2
		b1 <- n1 - k1 + a1
		b2 <- n2 - k2 + a2
		K <- beta(c1, b1)*beta(c2, b2)
		if(d >= (-1) & d <= 0){
			out <- try(beta(c1,b2)*F1(b2, c1 + c2 + b1 + b2 - 2, 1 - c2, c1 + b2, 1 + d, 1 - d^2, ...)*((-d)^(b1 + b2 - 1)*(1 + d)^(c1 + b2 - 1))/K, silent = TRUE)
		} else if(d > 0 & d <= 1){;
			out <- try(beta(c2, b1)*F1(b1, c1 + c2 + b1 + b2 - 2, 1 - c1, c2 + b1, 1 - d, 1 - d^2, ...)*(d^(b1 + b2 - 1)*(1 - d)^(c2 + b1 - 1))/K, silent = TRUE)
			} else out <- 0
		if(class(out) == "try-error") out <- 0
		out
	}
	TEMP <- Vectorize(ddiffprop2)
	temp <- TEMP(d = d, k1 = k1, k2 = k2, n1 = n1, n2 = n2, a1 = a1, a2 = a2,...)
    if(log) temp <- log(temp)
	temp
}

qdiffprop <- function(p, k1, k2, n1, n2, a1 = 0.5, a2 = 0.5, lower.tail = TRUE, log.p = FALSE,...)
{
    if ((a1 < 0 | a1 > 1) | (a2 < 0 | a2 > 1)) {
        warning("a1 and a2 should both be between 0 and 1 for this fiducial approach!", call. = TRUE)
    }
    if (log.p) p <- exp(p)
    if (lower.tail == FALSE) p <- 1 - p
	qdiffprop2 <- function(p, k1, k2, n1, n2, a1, a2,...)
	{
		if(p <= 0){
			out <- (-1)
		} else if(p >= 1){
			out <- 1
			} else{
				temp.fun <- function(d, p, k1, k2, n1, n2, a1, a2) p - integrate(Vectorize(ddiffprop), k1 = k1, k2 = k2, n1 = n1, n2 = n2, a1 = a1, a2 = a2, lower = -1, upper = d, ...)$value
				out <- try(uniroot(temp.fun, c(-1, 1), k1 = k1, k2 = k2, n1 = n1, n2 = n2, a1 = a1, a2 = a2, p = p)$root, silent = TRUE)
				if(class(out) == "try-error"){
					temp2 <- c(abs(temp.fun(-1,gamma,k1,k2,n1[1],n2[1],a1,a2)), abs(temp.fun(1,gamma,k1,k2,n1[1],n2[1],a1,a2)))
					out <- ifelse(which.min(temp2)==1,-1,1)
				}
			}
	out
	}
	TEMP <- Vectorize(qdiffprop2)
	temp <- TEMP(p, k1, k2, n1, n2, a1, a2, ...)
	temp
}

pdiffprop <- function(q, k1, k2, n1, n2, a1 = 0.5, a2 = 0.5, lower.tail = TRUE, log.p = FALSE,...)
{
    if ((a1 < 0 | a1 > 1) | (a2 < 0 | a2 > 1)) {
        warning("a1 and a2 should both be between 0 and 1 for this fiducial approach!", call. = TRUE)
    }
	x <- q
	pdiffprop2 <- function(x, k1, k2, n1, n2, a1, a2,...)
	{
		if(x <= (-1)){
			out <- 0
		} else if(x >= 1){
			out <- 1
			} else{
				out <- integrate(Vectorize(ddiffprop), k1 = k1, k2 = k2, n1 = n1, n2 = n2, a1, a2, lower = -1, upper = x, ...)$value
			}
				out
			}
	TEMP <- Vectorize(pdiffprop2)
	temp <- TEMP(x = x, k1 = k1, k2 = k2, n1 = n1, n2 = n2, a1, a2,...)
    if (lower.tail == FALSE) temp <- 1 - temp
    if (log.p) temp <- log(temp)
	temp
}

rdiffprop <- function(n, k1, k2, n1, n2, a1 = 0.5, a2 = 0.5)
{
    if ((a1 < 0 | a1 > 1) | (a2 < 0 | a2 > 1)) {
        warning("a1 and a2 should both be between 0 and 1 for this fiducial approach!", call. = TRUE)
    }
	out=rbeta(n, k1 + a1, n1 - k1 + a1) - rbeta(n, k2 + a2, n2 - k2 + a2)
	out
}
