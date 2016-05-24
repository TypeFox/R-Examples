# FREGAT (c) 2016 Gulnara R. Svishcheva, ICG SB RAS

NullMixedModel <- function(y, X, R, opt.method = 'optimize', ih2 = 0.3, eps = 1.e-04, H2est = TRUE, force = FALSE, ...) {
# dots to prevent 'unused argument' error

	#### y: phenotypes (nx1)
	#### X: covariate matrix (nxc)
	#### R: relationship matrix (nxn) (NOT kinship matrix)
	#### opt.method: three types of optimization ('optim', 'optimize', 'nlminb')
	#### ih2: starting value for h2
	#### H2est: logical value indicating whether 'h2' should be estimated.

	opt.method <- match.arg(opt.method, c('optimize', 'optim', 'nlminb'))

	#### likelihood function with regard to h2

	f_h2 <- function(h2) {
		Lam <- 1./(1. + h2*L) # (nx1) vector of eigenvalues for V^-1
		AAA <- t(tTX) %*% diag(Lam) # (mxn) matrix (AAA=t(tTX) %*% diag(Lam))
		Y <- tTy - (tTX %*% solve((AAA %*% tTX), AAA)) %*% tTy # (nx1) vector
		logLH <- log(mean(Lam*(Y*Y))) - mean(log(Lam)) # to minimize
		return(logLH) 
	}

	#### derivative of likelihood function with regard to h2

	df_h2 <- function(h2) {
		Lam <- 1./(1. + h2*L) # (nx1) vector of eigenvalues for V^-1
		AAA <- t(tTX) %*% diag(Lam) # (mxn) matrix
		Y <- tTy - (tTX %*% solve((AAA %*% tTX), AAA)) %*% tTy # (nx1) vector
		LamYY <- Lam*(Y*Y) # vector
		dfh2 <- mean(Lam*LamYY) - mean(Lam)*mean(LamYY) # it should be equal to 0
		return(dfh2) 
	}

	#### calculation of sigma2 and alpha (parameters under null model)

	Alpha_Sigma <- function(esth2) {
		Lam <- 1./(1. + esth2*L)
		AAA <- t(tTX) %*% diag(Lam)
		alpha <- solve((AAA %*% tTX), AAA) %*% tTy # regression coefficients for covariates
		Y <- tTy - tTX %*% alpha
		sigma2 <- mean(Lam*(Y*Y))
		return(list(alpha = alpha, sigma2 = sigma2))
	}

	#------------------------------------------------------------------------

	nnn <- length(y)

	if (H2est == TRUE){
		#### decomposition of relationship matrix R (R = T %*% diag(Lambda) %*% t(T))
		decomR <- eigen(R, symmetric = TRUE)
		Lambda <- decomR$values # Lambda is a vector of eigenvalues of R
		L <- Lambda - 1
		tT <- t(decomR$vectors) # tT is matrix of eigenvectors of R
		tTy <- tT %*% y # y is phenotype vector 
		tTX <- tT %*% X # X is covariate matrix

		#### methods for optimization
		if (opt.method == 'optim') {
			 A <- optim(par = ih2, fn = f_h2, gr = df_h2, lower = eps, upper = 1. - eps, method = "Brent")
			 Apar <- A$par
		} else if (opt.method == 'optimize') {
			B <- optimize(f = f_h2, interval = c(eps, 1.-eps))
			Apar <- B$minimum
		} else if (opt.method == 'nlminb') {
			C <- nlminb(objective = f_h2, start = ih2, lower = eps, upper = 1. - eps)
			Apar <- C$par
		}

		pr <- Alpha_Sigma(Apar)
		for(i in 1:length(pr)) assign(names(pr)[i], pr[[i]])

		df <- df_h2(Apar)
		LH = -0.5 * nnn * (log(2 * pi) + log(sigma2) + mean(log(1 + L*Apar)) + 1)

		sigma2 <- sigma2 * nnn / (nnn - 1)

		### fixed effects stat.significance

		CholCovInv <- try(chol(chol2inv(chol(sigma2 * (Apar * R + (1. - Apar) * diag(nnn)))))) # ma = t(chol(ma)) %*% chol(ma)
		if (is(CholCovInv, "try-error")) {
			if (!force) {
				stop("Covariance matrix has bad properties... Check the kinship matrix or try 'force = TRUE'")
			} else {
				CholCovInv <- chol(ginv(sigma2 * (Apar * R + (1. - Apar) * diag(nnn))))
			}
		}
		pheno <- CholCovInv %*% y # the tranformed phenotypes
		predicts <- CholCovInv %*% X # the tranformed covariates
		fit <- suppressWarnings(glm(pheno ~ predicts - 1, family = "gaussian"))
		a <- summary(fit)$coefficients
		rownames(a) <- gsub('predicts', '', rownames(a))
		rownames(a)[1] <- '(Intercept)'
		a[,1] <- alpha

		return(list(h2 = Apar, total.var = sigma2, alpha = a, df = df, logLH = LH))

	} else {
		fit <- suppressWarnings(glm(y ~ X - 1, family = "gaussian"))
		a <- summary(fit)$coefficients
		rownames(a) <- colnames(X)
#		rownames(a)[1] <- '(Intercept)'
		as.matrix(a[,1]) -> alpha
#		alpha <- (chol2inv(chol(t(X)%*%X)) %*% t(X)) %*% y
		yXa <- y - (X %*% alpha)
		sigma2 <- mean(yXa * yXa)
		LH <- -0.5 * nnn * (log(2 * pi) + log(sigma2) + 1) # total log-likelihood
		sigma2 <- sigma2 * nnn / (nnn - 1)
		rownames(a)[1] <- '(Intercept)'
		return(list(h2 = 0.0, total.var = sigma2, alpha = a, df = 0.0, logLH = LH)) # total var - var ostatkov..
	}

}

# ginv function (c) Package MASS version 7.3-23

ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
## Generalized Inverse of a Matrix
	dnx <- dimnames(X)
	if(is.null(dnx)) dnx <- vector("list", 2)
	s <- svd(X)
	nz <- s$d > tol * s$d[1]

	structure(
		if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
		dimnames = dnx[2:1])
}
