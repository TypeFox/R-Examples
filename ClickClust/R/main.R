
click.EM <- function(X, y = NULL, K, eps = 1e-10, r = 100, iter = 5, min.beta = 1e-3, min.gamma = 1e-3, scale.const = 1.0){

	if (K < 1) stop("Wrong number of mixture components K...\n")
	if (eps <= 0) stop("Wrong value of eps...\n")
	if (r < 1) stop("Wrong number of random restarts r...\n")
	if (iter < 1) stop("Wrong number of iterations iter...\n")
	if (min.beta < 0) stop("Wrong lower bound min.beta...\n")
	if (min.gamma < 0) stop("Wrong lower bound min.gamma...\n")
	if (scale.const <= 0) stop("Wrong value of scale.const...\n")

	xdims <- dim(X)
	p <- xdims[1]
	n <- xdims[3]

	x1 <- as.vector(X)
	id <- rep(0, n)

	alpha <- rep(0, K)
	beta <- rep(0, K*p)
	gamma <- rep(-1, p*p*K)
	z <- rep(0, n*K)
	l <- c(0, 0)

	if (!is.null(y)){

		y <- y - 1

		Q <- .C("runClickClust_", p1 = as.integer(p), K1 = as.integer(K), n1 = as.integer(n), x1 = as.integer(x1), y = as.integer(y), alpha = as.double(alpha), beta1 = as.double(beta), Pi1 = as.double(gamma), gamma1 = as.double(z), id = as.integer(id), e1 = as.double(eps), r1 = as.integer(r), shortem1 = as.integer(iter), l = as.double(l), lowbeta1 = as.double(min.beta), lowPi1 = as.double(min.gamma), scaleconst1 = as.double(scale.const), PACKAGE = "ClickClust")

	} else {

		Q <- .C("runClickClust", p1 = as.integer(p), K1 = as.integer(K), n1 = as.integer(n), x1 = as.integer(x1), alpha = as.double(alpha), Pi1 = as.double(gamma), gamma1 = as.double(z), id = as.integer(id), e1 = as.double(eps), r1 = as.integer(r), shortem1 = as.integer(iter), l = as.double(l), lowPi1 = as.double(min.gamma), scaleconst1 = as.double(scale.const), PACKAGE = "ClickClust")

	}

	a <- array(Q$Pi1, c(K, p, p))
	b <- array(NA, c(p, p, K))
	for (i in 1:K){
		b[,,i] <- t(a[i,,])
	}
	b[b == -1] <- NA

	if (!is.null(y)){

		ret <- list(z = t(matrix(Q$gamma1, nrow = K)), id = Q$id + 1, alpha = Q$alpha, beta = matrix(Q$beta, ncol = p, byrow = TRUE), gamma = b, logl = Q$l[1], BIC = Q$l[2])

	} else {

		ret <- list(z = t(matrix(Q$gamma1, nrow = K)), id = Q$id + 1, alpha = Q$alpha, gamma = b, logl = Q$l[1], BIC = Q$l[2])

	}

	class(ret) <- "EM"
	return(ret)

}


print.EM <- function(x, ...){
        K <- length(x$alpha)
        p <- dim(x$gamma)[1]
	cat("\nK = ", K,
	    ", p = ", p,
            ", logl = ", x$logl,
            ", BIC = ", x$BIC,
            "\n", sep = "")
        cat("\nCluster sizes:")
        print(table(x$id))
        cat("\nid: \n")
        print(x$id)
        cat("\nalpha: \n")
        print(x$alpha)
	if (!is.null(x$beta)){
		cat("\nbeta: \n")
        	print(x$beta)
	}
        cat("\nUse $ to access:\n\t-transition probability matrices (gamma)\n\t-posterior probabilities (z)\n")
	invisible()
} # End of print.EM().


summary.EM <- function(object, ...){
        K <- length(object$alpha)
        p <- dim(object$gamma)[1]
	cat("\nK = ", K,
	    ", p = ", p,
            ", logl = ", object$logl,
            ", BIC = ", object$BIC,
            "\n", sep = "")
        cat("\nCluster sizes:")
        print(table(object$id))
	invisible()
} # End of summary.EM().


print.search <- function(x, ...){
        K <- length(x$alpha)
        p <- length(x$states)
        d <- max(x$states)
	cat("\nK = ", K,
	    ", p = ", p,
	    ", d = ", d,
            ", logl = ", x$logl, sep = "")
	if (is.null(x$BIC)) cat(", AIC = ", x$AIC, "\n", sep = "")
		else cat(", BIC = ", x$BIC, "\n", sep = "")
        cat("\nStates:")
        print(x$states)
        cat("\nCluster sizes:")
        print(table(x$id))
        cat("\nid: \n")
        print(x$id)
        cat("\nalpha: \n")
        print(x$alpha)
        cat("\nUse $ to access:\n\t-transition probability matrices (gamma)\n\t-posterior probabilities (z)\n")
	invisible()
} # End of print.search().


summary.search <- function(object, ...){
        K <- length(object$alpha)
        p <- length(object$states)
        d <- max(object$states)
	cat("\nK = ", K,
	    ", p = ", p,
	    ", d = ", d,
            ", logl = ", object$logl, sep = "")
	if (is.null(object$BIC)) cat(", AIC = ", object$AIC, "\n", sep = "")
		else cat(", BIC = ", object$BIC, "\n", sep = "")
        cat("\nStates:")
        print(object$states)
        cat("\nCluster sizes:")
        print(table(object$id))
	invisible()
} # End of summary.search().







