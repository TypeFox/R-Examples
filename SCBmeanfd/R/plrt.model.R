plrt.model <- function(x, y, model, verbose = FALSE)
{

	# Check arguments
	spacing <- diff(x)
	if (any(spacing < 0) || !isTRUE(all.equal(min(spacing), max(spacing))))
		stop("x must be a uniform grid") 
	p  <- length(x)
	if (ncol(y) != p) 
		stop("length(x) and ncol(y) must be equal")
	PolyModel <- (length(model) == 1)
	if (!PolyModel && nrow(model) != p) {
 		stop("'model' must be an integer or a matrix such that nrow(model) == length(x)")
 	}
 	
	# Data mean
	ymean <- colMeans(y)
	n 	  <- nrow(y) 

	# Parametric fit	
	fit0  <- if (PolyModel) { lm(ymean ~ poly(x, degree = model))
			 } else lm(ymean ~ model + 0)
		
	# Projection matrix onto parametric residuals	
	I  <- diag(p)
	X  <- model.matrix(fit0)
	M0 <- I - X %*% solve(crossprod(X)) %*% t(X)
	e <- residuals(fit0)	# = M0 %*% ymean
	
	# Projection matrix onto nonparametric residuals
	h  <- 2 * (max(x) - min(x)) / p	
	f  <- function(y) locpoly(x = x, y = y, bandwidth = h, gridsize = p)$y
	S  <- apply(I, 2, f)
	M1 <- crossprod(I - S)

	# Test statistic
	rss0  <- as.numeric(crossprod(e))
	rss1  <- as.numeric(crossprod(e, M1 %*% e))
	Ftest <- rss0 / rss1 - 1
	
	# Chi square approximation to distribution of pseudo-likelihood ratio
	V 	  <- cov(as.matrix(y)) / n
	A 	  <- V %*% M0 %*% (I - (1 + Ftest) * M1) %*% M0	
	A2    <- A %*% A
	A3    <- A2 %*% A
	k1 	  <- sum(diag(A))
	k2    <- 2 * sum(diag(A2))
	k3 	  <- 8 * sum(diag(A3))
	a 	  <- k3 / (4 * k2)
	b 	  <- 8 * k2^3 / k3^2
	c 	  <- k1 - a * b
	
	# p-value
	pval  <- pchisq(-c/a, b, lower.tail = (a < 0))
		
	if (verbose) {
		cat("\nPseudo-Likelihood Ratio Test\n")
		cat ("Model for the mean function: ")
		if (length(model) == 1) {
			if (model == 0) { cat ("zero\n") 
			} else if (model == 1) { cat("linear\n")
			} else cat("polynomial of degree <=", model,"\n")
		} 
		else cat("function space of dimension", ncol(model),"\n")	
		cat("Bandwidth:", round(h, 4), "\nTest statistic and p value\n")
		df <- data.frame(F = Ftest, p = pval)
		print(df, digits = 4, row.names = FALSE)
		return(invisible(pval))
	}
	return(pval)
}
