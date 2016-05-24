bigRR.default <-
function (formula = NULL, y, X, Z, data = NULL, shrink = NULL, weight = NULL,
    family = gaussian(link = identity), lambda = NULL, impute = FALSE, tol.err = 1e-6, 
	tol.conv = 1e-8, only.estimates = FALSE, GPU = FALSE, ...) 
{

    Call <- match.call()
    if (!(is.matrix(X))) 
        stop("X should be a matrix.")
    if (!(is.matrix(Z))) 
        stop("Z should be a matrix.")
    if (!(is.vector(y))) 
        stop("y should be a vector.")
	#if (GPU & !require(gputools, quietly = TRUE)) stop('Package gputools is needed for using GPU.\n')
	if (GPU) stop('GPU option is only available in R-Forge versions of "bigRR".\n')
	if (any(is.na(y))) {
		naidx <- which(is.na(y))
		y <- y[-naidx]
		X <- X[-naidx,]
		Z <- Z[-naidx,]
	}
	if (impute) {
		if (any(is.na(Z))) {
			cat('Imputing missing values...')
			nacolidx <- which(is.na(colSums(Z)))
			for (j in nacolidx) {
				naidx <- which(is.na(Z[,j]))
				Z[naidx,j] <- sample(Z[-naidx,j], length(naidx), TRUE)
			}
			cat('Done.\n')
		} else {
			cat("NOTE: no missing value exists, no need to impute.\n")
		}
	}
    N <- n <- nrow(X)
    p <- ncol(X)
    k <- ncol(Z)
    if (N != nrow(Z) | N != length(y)) 
        stop("Sizes of y, X, and Z are not all equal.")
    if (is.null(weight)) w <- rep(1, k) else w <- weight

    #G <- crossprod(sqrt(w)*t(Z)) ## bug fixed 111201 -- Xia
	wZt <- sqrt(w)*t(Z)
	if (!GPU) G <- crossprod(wZt) #else G <- gpuMatMult(t(wZt), wZt)

	############ Bending to allow for p<n problems -- Lars (Xia added SVD)
	if (k < n) {
		eigen.values <- eigen(G)$values
		min.eigen <- min(eigen.values)
		if (min.eigen < tol.err) G <- G + diag(N)*(abs(min.eigen) + tol.err) 
	}
	############
	#invG <- solve(G)
	#L <- t(chol(G))
	svdG <- svd(G)
	L <- svdG$u %*% diag(sqrt(svdG$d))
	invG <- tcrossprod(svdG$v %*% diag(1/svdG$d), svdG$u)
    phi0 <- sa0 <- 1
    if (is.null(lambda)) {
        hm <- hglm(y = y, X = X, Z = L, family = family, conv = tol.conv) ## checked with old emme code, conv = 1e-6 removed -- Xia
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
	if (!GPU) tZinvG <- crossprod(Z, invG) #else tZinvG <- gpuMatMult(t(Z), invG)
    u <- (w*tZinvG)%*%a
    qu <- GCV <- NULL
    if (!only.estimates) {
        C <- rbind(cbind(crossprod(X, X), crossprod(X, L)), cbind(crossprod(L, X), G + diag(N)*phi/sa))
        C22 <- solve(C)[(p + 1):(p + N), (p + 1):(p + N)]*phi
		if (!GPU) transf <- (w*tZinvG) %*% L #else transf <- gpuMatMult((w*tZinvG), L)
        qu <- hat.transf(C22, transf, vc = sa, w, k, N, tol.err = tol.err, GPU = GPU)
		adjsmahat <- adjbighat <- 0
		if (any(qu < tol.err) | any (qu > 1 - tol.err)) {
			adjsmahat <- sum(qu < tol.err)
			adjbighat <- sum(qu > (1 - tol.err))
        	qu[qu < tol.err] <- tol.err
        	qu[qu > (1 - tol.err)] <- 1 - tol.err
		}
	    if (family$family == "gaussian") GCV <- sum(hm$resid^2)/((n - sum(hm$hv[1:n]))^2)
    }
    result <- list(phi = phi, lambda = hm$varRanef, beta = hm$fixef, hglm = hm,
        u = u, leverage = qu, GCV = GCV, Call = Call, y = y, X = X, Nsmallhat = adjsmahat, Nbighat = adjbighat)
    class(result) <- "bigRR"
    return(result)
}

