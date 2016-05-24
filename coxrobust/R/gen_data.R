
gen_data <- function(n, beta, cont = 0, p.censor = 0) {

	beta <- as.double(beta)
	m <- length(beta)

    z <- array(rnorm(n*m), c(n,m))

    time <- rexp(n) / exp( z %*% beta )

	status <- sample(c(0,1), n, replace = TRUE, prob = c(p.censor, 1-p.censor))

	ncont <- floor(m*n*cont)
	if ( ncont > 0 ) {
		z <- as.double(t(z))
    	z[1:ncont] <- 2*rnorm(ncont) + 1
		z <- matrix(z, n, m, byrow = TRUE)
	}

    if ( ncol(z) == 1 ) {
        gdata <- data.frame(time, status, X1 = z)
    } else {
        gdata <- data.frame(time, status, z)
    }
    
    return(gdata)

}
