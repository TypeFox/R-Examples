Manly.EM <- function(X, id = NULL, la = NULL, tau = NULL, Mu = NULL, S = NULL, tol = 1e-5, max.iter = 1000){


	if (tol <= 0) stop("Wrong value of tol...\n")
	if (max.iter < 1) stop("Wrong number of iterations iter...\n")


	xdims <- dim(X)
	n <- xdims[1]
	p <- xdims[2]


	if(n < 1) stop("Wrong number of observations n...\n")
	if(p < 1) stop("Wrong dimensionality p...\n")



	if(is.null(id) && (is.null(Mu) || is.null(tau) || is.null(S))) stop("Must provide one initialization method...\n")
	if(!is.null(id)){


		K <- max(id)	

		if(K < 1) stop("Wrong number of mixture components K...\n")
		if(is.null(la)){
			la <- matrix(0.0, K, p)
		}
		if(K != dim(la)[1]) stop("Inconsistent number of mixture components K...\n")	
		if(p != dim(la)[2]) stop("Inconsistent dimensionality p...\n")


		x1 <- as.vector(t(X))
		la1 <- as.vector(t(la))
		gamma1 <- rep(0, n*K)
		ll <- rep(0, 3)
		misc_int <- c(p, n, K, max.iter)
		misc_double <- c(tol, 0.0, 0.0)
		conv <- rep(0, 2)
		tau <- rep(0, K)
		Mu1 <- rep(0, K*p)
		S1 <- rep(0, K*p*p)


		result <- .C("run_Manly", x1 = as.double(x1), id = as.integer(id), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), la1 = as.double(la1), tau = as.double(tau), Mu1 = as.double(Mu1), S1 = as.double(S1), gamma1 = as.double(gamma1), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "ManlyMix")
	}
	else{
		K <- length(tau)
		if(is.null(la)){
			la <- matrix(0.0, K, p)
		}

		equal.K <- c(dim(la)[1], dim(Mu)[1], dim(S)[3])
		equal.p <- c(dim(la)[2], dim(Mu)[2], dim(S)[1], dim(S)[2])

		if (K < 1) stop("Wrong number of mixture components K...\n")
		if ((K != equal.K[1]) || (K != equal.K[2]) || (K != equal.K[3])) stop("Inconsistent number of mixture components K...\n")
		if ((p != equal.p[1]) || (p != equal.p[2]) || (p != equal.p[3]) || (p != equal.p[4])) stop("Inconsistent number of dimensionality p...\n")


		x1 <- as.vector(t(X))
		la1 <- as.vector(t(la))
		gamma1 <- rep(0, n*K)
		ll <- rep(0, 3)
		misc_int <- c(p, n, K, max.iter)
		misc_double <- c(tol, 0.0, 0.0)
		conv <- rep(0, 2)
		id <- rep(0, n)
		Mu1 <- as.vector(t(Mu))
		S1 <- as.vector(S)


		result <- .C("run_Manly2", x1 = as.double(x1), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), la1 = as.double(la1), tau = as.double(tau), Mu1 = as.double(Mu1), S1 = as.double(S1), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv),PACKAGE = "ManlyMix")

	}

	ind <- result$la1==0
	
	M <- K - 1 + 2 * K * p + K * p * (p + 1) / 2 - sum(ind)
	BIC <- Manly.bic(result$ll[1], n, M)

	
	ret <- list(la = matrix(result$la1, nrow=K, byrow=TRUE), tau = result$tau, Mu = matrix(result$Mu1, nrow = K, byrow = TRUE), S = array(result$S1, dim = c(p, p, K)), gamma = matrix(result$gamma1, nrow = n, byrow = TRUE), id = result$id, ll = result$ll[1], bic = BIC, iter = result$conv[1], flag = result$conv[2])

	class(ret) <- "ManlyMix"
	if(result$conv[2] == 1) {
		warning("The EM algorithm does not converge...\n")
		ret <- list(la = NULL, tau = NULL, Mu = NULL, S = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = result$conv[1], flag = result$conv[2])
	}	
	return(ret)



}


Manly.bic <- function(logl, n, M){
	return(-2 * logl + M * log(n))
}





