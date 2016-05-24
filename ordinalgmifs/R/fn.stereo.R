fn.stereo <-
function(par, w, x, beta, y, k, Ymat) {
	if (dim(w)[2]!=0) {
	    theta<- as.vector(par[(k-1+k-2+1):length(par)]) 
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
    eta<-matrix(0, ncol=k-1, nrow=length(y))
    eta[,1] <- exp(par[1] + Xb)
	for (i in 2:(k-1)) {
		eta[,i]<- exp(par[i] + par[i+k-2]*Xb)
	}
	eta.sum <- .Call("do_matrix_sum_rows", eta)
	pik <- 1 - eta.sum/(1+eta.sum) 
	pi <- matrix(0,ncol=k,nrow=length(y))
	pi[,k] <- pik
	pi[,1:(k-1)] <- eta*pik
	mult<-Ymat*log(pi)
	mult.vec <- .Call("do_matrix_sum_rows", mult)
	loglik <- .Call("do_vector_sum", mult.vec)
	-loglik
}
