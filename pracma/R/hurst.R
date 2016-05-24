##
##  h u r s t . R  Hurst Exponent
##


.hurstrs <- function(x) {
	# half intervals of indices
	half <- function(N) sort(c(N, N[-length(N)]+((diff(N)+1)%/%2)))
	# define the R/S scale
	rscalc <- function(x) {
		n <- length(x); y <- cumsum(x - mean(x))
		R <- diff(range(y)); S <- sd(x)
		return(R/S)
	}
	# set initial values
	X <- c(length(x))
	Y <- c(rscalc(x))
	N <- c(0, length(x) %/% 2, length(x))
	# compute averaged R/S for halved intervals
	while ( min(diff(N)) >= 8 ) {
		xl <- c(); yl <- c()
		for (i in 2:length(N)) {
			rs <- rscalc(x[(N[i-1]+1):N[i]])
			xl <- c(xl, N[i]-N[i-1])
			yl <- c(yl, rs)
		}
		X <- c(X, mean(xl))
		Y <- c(Y, mean(yl))
		# next step
		N <- half(N)
	}
	# apply linear regression
	rs_lm <- lm(log(Y) ~ log(X))
	return(unname(coefficients(rs_lm)[2]))
}


hurstexp <- function(x, d = 50, display = TRUE) {
    stopifnot(is.numeric(x), is.numeric(d))
    d <- max(2, floor(d[1]))
    N <- length(x)
    if (N %% 2 != 0) {
        x <- c(x, (x[N-1] + x[N])/2)
        N <- N + 1
    }

    # Calculate simple R/S
    rssimple <- function(x){
        n <- length(x)
        y <- x - mean(x)
        s <- cumsum(y)
        rs <- (max(s) - min(s)) / sd(x)
        log(rs) / log(n)
    }

    # Calculate empirical R/S
    rscalc <- function(z, n) {
        m <- length(z)/n
        y <- matrix(x, n, m)
        e <- apply(y, 2, mean)
        s <- apply(y, 2, std)
        for (i in 1:m) y[, i] <- y[, i] - e[i]
        y <- apply(y, 2, cumsum)
        mm <- apply(y, 2, max) - apply(y, 2, min)
        return( mean(mm/s) )
    }
    divisors <- function(n, n0 = 2) {
        n0n <- n0:floor(n/2)
        dvs <- n0n[n %% n0n == 0]
        return(dvs)
    }

    # Find the optimal vector d
    N <- length(x); dmin <- d
    N0 <- min(floor(0.99 * N), N-1)
    N1 <- N0; dv <- divisors(N1, dmin)
    for (i in (N0+1):N) {
        dw <- divisors(i, dmin)
        if (length(dw) > length(dv))
            N1 <- i; dv <- dw
    }
    OptN <- N1; d <- dv
    x <- x[1:OptN]

    N <- length(d)
    RSe <- ERS <- numeric(N)
    for (i in 1:N)
        RSe[i] <- rscalc(x, d[i])

    # Compute corrected theoretical E(R/S)
    for (i in 1:N) {
        n <- d[i]
        K <- c((n-1):1)/c(1:(n-1))
        ratio <- (n-0.5)/n * sum(sqrt(K))
        if (n > 340) ERS[i] <- ratio/sqrt(0.5*pi*n)
        else         ERS[i] <- (gamma(0.5*(n-1))*ratio) / (gamma(0.5*n)*sqrt(pi))
    }

    # Compute the Hurst exponent as the slope on a loglog scale
    ERSal <- sqrt(0.5*pi*d)
    Pal <- polyfit(log10(d), log10(RSe - ERS + ERSal), 1)
    Hal <- Pal[1]

    # Calculate the empirical and theoretical Hurst exponents
    Pe <- polyfit(log10(d), log10(RSe), 1)
    He <- Pe[1]
    P  <- polyfit(log10(d), log10(ERS), 1)
    Ht <- P[1]

    Hs   <- rssimple(x)
    Hrs <- .hurstrs(x)

    if (display) {
        cat("Simple R/S Hurst estimation:        ", Hs,  "\n")
        cat("Corrected R over S Hurst exponent:  ", Hrs, "\n")
        cat("Empirical Hurst exponent:           ", He,  "\n")
        cat("Corrected empirical Hurst exponent: ", Hal, "\n")
        cat("Theoretical Hurst exponent:         ", Ht,  "\n")
        invisible(list(Hs = Hs, Hrs = Hrs, He = He, Hal = Hal, Ht = Ht))
    } else {
        return(list(Hs = Hs, Hrs = Hrs, He = He, Hal = Hal, Ht = Ht))
    }
}
