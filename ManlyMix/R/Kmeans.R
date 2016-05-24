

Manly.Kmeans <- function(X, id = NULL, la = NULL, Mu = NULL, S = NULL, tol = 1e-5, max.iter = 1000){


	if (tol <= 0) stop("Wrong value of tol...\n")
	if (max.iter < 1) stop("Wrong number of iterations iter...\n")


	xdims <- dim(X)
	n <- xdims[1]
	p <- xdims[2]


	if(n < 1) stop("Wrong number of observations n...\n")
	if(p < 1) stop("Wrong dimensionality p...\n")



	if(is.null(id) && (is.null(Mu) || is.null(S))) stop("Must provide one initialization method...\n")
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

		ll <- rep(0, 3)
		misc_int <- c(p, n, K, max.iter)
		misc_double <- c(tol, 0.0, 0.0)
		conv <- rep(0, 2)

		Mu1 <- rep(0, K*p)
		S1 <- rep(0, K)


		result <- .C("run_Manlyk", x1 = as.double(x1), id = as.integer(id), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), la1 = as.double(la1), Mu1 = as.double(Mu1), S1 = as.double(S1), conv = as.integer(conv), PACKAGE = "ManlyMix")
	}
	else{
		K <- dim(Mu)[1]
		if(is.null(la)){
			la <- matrix(0.0, K, p)
		}

		equal.K <- c(dim(Mu)[1], length(S))
		equal.p <- c(dim(la)[2], dim(Mu)[2])

		if (K < 1) stop("Wrong number of mixture components K...\n")
		if ((K != equal.K[1]) || (K != equal.K[2])) stop("Inconsistent number of mixture components K...\n")
		if ((p != equal.p[1]) || (p != equal.p[2])) stop("Inconsistent number of dimensionaltiy p...\n")


		x1 <- as.vector(t(X))
		la1 <- as.vector(t(la))
	
		ll <- rep(0, 3)
		misc_int <- c(p, n, K, max.iter)
		misc_double <- c(tol, 0.0, 0.0)
		conv <- rep(0, 2)
		id <- rep(0, n)
		Mu1 <- as.vector(t(Mu))
		S1 <- as.vector(S)


		result <- .C("run_Manlyk2", x1 = as.double(x1), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), la1 = as.double(la1), Mu1 = as.double(Mu1), S1 = as.double(S1), id = as.integer(id), conv = as.integer(conv),PACKAGE = "ManlyMix")

	}
	

	ret <- list(la = matrix(result$la1, nrow=K, byrow=TRUE), Mu = matrix(result$Mu1, nrow = K, byrow = TRUE), S = result$S1, id = result$id, iter = result$conv[1], flag = result$conv[2])
	if(result$conv[2] == 1) {
		warning("The CEM algorithm does not converge...\n")
		ret <- list(la = NULL, Mu = NULL, S = NULL, id = NULL, iter = result$conv[1], flag = result$conv[2])
	}	
	class(ret) <- "ManlyMix"
	
	
	return(ret)



}









