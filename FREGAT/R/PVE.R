# FREGAT (c) 2016 Gulnara R. Svishcheva & Nadezhda M. Belonogova, ICG SB RAS

tau.fun <- function(a, lam, zzz) { # reml estimate

	tau2 <- a / (1 - a)
	sum(log(1 + tau2 * lam)) - tau2 * sum((zzz ^ 2)/(1 + tau2 * lam))

}

tau3 <- function(param, X, R, Ker, y, I) { # full likelihood (3 case)
	sigma2 <- param[1]                # total variance
	hk2 <- param[2]                   # total heritability h2
	dd <- param[3]                    # hk2*dd :  hk2 * dd - proportion of variance by relationship (R)
	X <- as.matrix(X)                 # matrix of covariates
	alpha <- param[4:(3 + dim(X)[2])] # intercept and alphas
	ya <- y - X %*% alpha             # centralized phenotypes
	kkk <- dim(R)[1] / sum(diag(Ker))   # coefficient for normalization of Ker
	Cov <- (sigma2 * hk2 * dd) * R + sigma2 * hk2 * (1. - dd) * kkk * Ker + (sigma2 * (1. - hk2)) * I
	eig <- eigen(Cov, symmetric = TRUE)
	neg.Log2.LH <- sum((log(eig$value[eig$value > 0]))) + sum(ya * (ginv(Cov) %*% ya))
	# to minimize
	return(neg.Log2.LH)
}

tau3.old <- function(param, X, R, Ker, y, I) { # full likelihood

	sigma2 <- param[1] # total variance
	hk2 <- param[2]
	dd <- param[3]
	X <- as.matrix(X)
	alpha <- param[4:(3 + dim(X)[2])]
	ya <- y - X %*% alpha
	env2 <- 1 - hk2 # proportion of variance by environment
	Cor <- (hk2 * dd) * R + (hk2 * (1. - dd)) * Ker + env2 * I
	# hk2 = (sigmah2 + tau2) / sigma2
	# dd = sigmah2 / (sigmah2 + sigmaK2)
	Cov <- sigma2 * Cor
	eig <- eigen(Cov, symmetric = TRUE)
	neg.Log2.LH <- sum((log(eig$value[eig$value > 0]))) + sum(ya * (ginv(Cov) %*% ya)) # to minimize
	# ginv(Cov) %*% ya ~= solve(Cov, ya)
	return(neg.Log2.LH)

}
