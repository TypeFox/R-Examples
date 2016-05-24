fn.acat <-
function(par, w, x, beta, y, k, Ymat) {
	if (dim(w)[2]!=0) {
    	theta <- par[k:length(par)]
		if (is.null(x)) {
    		Xb <- w%*%theta
    	} else if (!is.null(x)) {
    		Xb <- w%*%theta + x%*%beta
    	}
    } else if (!is.null(x)) {
    		Xb <- x%*%beta
    } else {
    		Xb <-0
    }
	eta <- matrix(0, ncol=k-1, nrow=dim(w)[1])
	for (i in 1:(k-1)) {
		eta[,i] <- par[i] + Xb
	}
	ceta <- .Call("do_row_cumsum", eta)
	eta.cumsum <- matrix(ceta,
		nrow=nrow(eta), byrow=T)
	numer <- rep(0, dim(eta.cumsum)[1])
	for (i in 1:dim(eta.cumsum)[2]) {
		numer <- numer + exp(eta.cumsum[,i])
	}
	pi <- matrix(0, ncol=k, nrow=dim(w)[1])
	pi[,1] <- 1 - numer / (1 + numer)
	for (i in 2:k) {
		pi[,i] <- exp(eta.cumsum[,i-1] + log(pi[,1]))
	}
	#loglik <- sum(apply(Ymat * log(pi), 1, sum))
	loglik <- sum(.Call("do_row_product_sums", Ymat, log(pi)))
	-loglik
}
