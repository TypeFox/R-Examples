mlfact <- function (Sigma, nfact){    
	# Phil Chalmers, December 7, 2011 
	# Following code substantially modified from stats::factanal function.
	FAfn <- function(Psi, S, q) {
		sc <- diag(1/sqrt(Psi))
		Sstar <- sc %*% S %*% sc
		E <- eigen(Sstar, symmetric = TRUE, only.values = TRUE)
		e <- E$values[-(1L:q)]
		e <- sum(log(e) - e) - q + nrow(S)
		-e
	}
	FAgr <- function(Psi, S, q) {
		sc <- diag(1/sqrt(Psi))
		Sstar <- sc %*% S %*% sc
		E <- eigen(Sstar, symmetric = TRUE)
		L <- E$vectors[, 1L:q, drop = FALSE]
		load <- L %*% diag(sqrt(pmax(E$values[1L:q] - 1, 0)), 
			q)
		load <- diag(sqrt(Psi)) %*% load
		g <- load %*% t(load) + diag(Psi) - S
		diag(g)/Psi^2
	}
	p <- ncol(Sigma)
	pars <- (1 - 0.5 * nfact/p)/diag(solve(Sigma))	
	res <- optim(pars, FAfn, FAgr, method = "L-BFGS-B", lower = .005, 
		upper = 1, control = c(list(fnscale = 1, parscale = rep(0.01, 
			length(pars)))), q = nfact, S = Sigma, hessian = TRUE)
	sc <- diag(1/sqrt(res$par))
	Sstar <- sc %*% Sigma %*% sc
	E <- eigen(Sstar, symmetric = TRUE)
	L <- E$vectors[, 1L:nfact, drop = FALSE]
    load <- L %*% diag(sqrt(pmax(E$values[1L:nfact] - 1, 0)), nfact)
    res$loadings <- diag(sqrt(res$par)) %*% load
	res
}
