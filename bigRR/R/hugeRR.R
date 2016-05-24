hugeRR <-
function(y, X, Z.name, Z.index, weight = NULL, family = gaussian(link = identity), 
		 lambda = NULL, only.estimates = FALSE, tol.err = 1e-6, tol.conv = 1e-8, save.cache = FALSE, ...) {
	Call <- match.call()
	if (!(is.matrix(X))) 
		stop("X should be a matrix.")
	if (!(is.vector(y))) 
		stop("y should be a vector.")
	N <- n <- nrow(X)
	p <- ncol(X)
	k <- 0
	if (!file.exists(paste('G_stuff_', sum(weight), '.RData', sep = ''))) {
		cat('Creating kinship matrix ...\n')
		for (i in Z.index) {
			subZ <- databel(paste(Z.name, i, sep = ''))
			subZ <- databel2matrix(subZ, 1:nrow(subZ), 1:ncol(subZ))
			if (i == Z.index[1]) G <- matrix(0, nrow(subZ), nrow(subZ))
			if (!is.null(weight)) subZ <- t(sqrt(weight[(k + 1):(k + ncol(subZ))])*t(subZ))
			G <- G + tcrossprod(subZ)
			k <- k + ncol(subZ)
			cat(i, '\t')
		}
		cat('\n')
		if (N != nrow(G) | N != length(y)) 
			stop("Sizes of y, X, and Z are not all equal.")
	
		cat('Cholesky-decomposition & inversion ...\n')
		svdG <- svd(G)
		L <- svdG$u %*% diag(sqrt(svdG$d))
		invG <- tcrossprod(svdG$v %*% diag(1/svdG$d), svdG$u)
		if (save.cache) save(G, L, invG, file = paste('G_stuff_', sum(weight), '.RData', sep = ''))
	} else {
		load(paste('G_stuff_', sum(weight), '.RData', sep = ''))
	}
	
	cat('Fitting HGLM ...\n')
	if (is.null(lambda)) {
		hm <- hglm(y = y, X = X, Z = L, family = family, conv = tol.conv) 
	}
	else {
		start.beta = c(rep(0, p))
		start.v = c(rep(0.01, n))
		start.lambda = lambda
		start.sigma2e = 1
		cat("Only 1 iteration applied for fixed lambda")
		hm <- hglm(y = y, X = X, Z = L, family = family, startval = c(start.beta, 
						start.v, start.lambda, start.sigma2e), maxit = 1)
	}
	phi <- as.numeric(hm$varFix)
	sa <- as.numeric(hm$varRanef)
	a <- L%*%hm$ranef
	#tZinvG <- crossprod(Z, invG)
	k2 <- 0
	u <- GCV <- qu <- NULL
	cat('Calculating effects and leverages ...\n')
	for (i in Z.index) {
		if (!file.exists(paste('tZinvG', i, '_', sum(weight), '.RData', sep = ''))) {
			subZ <- databel(paste(Z.name, i, sep = ''))
			subZ <- databel2matrix(subZ, 1:nrow(subZ), 1:ncol(subZ))
			tZinvG <- crossprod(subZ, invG)
			if (!is.null(weight)) tZinvG <- weight[(k2 + 1):(k2 + ncol(subZ))]*tZinvG
			if (save.cache) save(tZinvG, file = paste('tZinvG', i, '_', sum(weight), '.RData', sep = ''))
		} else {
			load(paste('tZinvG', i, '_', sum(weight), '.RData', sep = ''))
		}
		u <- c(u, tZinvG %*% a)
		if (!only.estimates) {
			C <- rbind(cbind(crossprod(X, X), crossprod(X, L)), cbind(crossprod(L, X), G + diag(N)*phi/sa))
			C22 <- solve(C)[(p + 1):(p + N), (p + 1):(p + N)]*phi
			transf <- tZinvG %*% L			
			#qu <- hat.transf(C22, transf, vc = sa, w, k, N, tol.err = tol.err, GPU = GPU)
			middle <- diag(N) - C22/sa
			svdmid <- svd(middle)
			M <- t(svdmid$u %*% diag(sqrt(svdmid$d)))
			if (!is.null(weight)) {
				qu <- c(qu, 1 - colSums(tcrossprod(M, transf)**2)/weight[(k2 + 1):(k2 + nrow(tZinvG))])
			} else {
				qu <- c(qu, 1 - colSums(tcrossprod(M, transf)**2))
			}
		}
		k2 <- k2 + nrow(tZinvG)
		cat(i, '\t')
	}
	cat('\nDone.\n')
	qu[qu < tol.err] <- tol.err
	qu[qu > (1 - tol.err)] <- 1 - tol.err
	if (!only.estimates & family$family == "gaussian") GCV <- sum(hm$resid^2)/((n - sum(hm$hv[1:n]))^2)
	
	result <- list(phi = phi, lambda = hm$varRanef, beta = hm$fixef, hglm = hm,
			u = u, leverage = qu, GCV = GCV, Call = Call, y = y, X = X)
	class(result) <- "bigRR"
	return(result)
}

