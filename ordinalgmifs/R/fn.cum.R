fn.cum <-
function (theta, w, x, beta, y, k, levels, Ymat, link) {
    alpha <- theta[1:(k-1)]
	if (dim(w)[2]!=0) {
    	zeta <- theta[k:length(theta)]
		if (is.null(x)) {
    		Xb <- w%*%zeta
    	} else if (!is.null(x)) {
    		Xb <- w%*%zeta + x%*%beta
    	}
    } else if (!is.null(x)) {
    		Xb <- x%*%beta
    } else {
    		Xb <-0
    }
    if (link=="logit") {
		z <- matrix(ncol=k-1, nrow=dim(w)[1])
    	for (i in 1:(k-1)) {
	    	z[,i] <- alpha[i] + Xb
		}
  		pi.z <- matrix(ncol=k, nrow=dim(w)[1])
		pi.z <- .Call("do_exp", k, z, pi.z)
#		pi <- .Call("do_row_product_sums", pi.z, Ymat)
#		loglik <- .Call("do_vector_sum", log(pi))
		pi <- pi.z*Ymat		
		pi <- .Call("do_row_products", pi)
		loglik <- .Call("do_vector_sum", log(pi))	
#		loglik <- .Call("do_vector_sum", log(pi.z)*Ymat)
	} else if (link=="probit") {
	    G.mat <- matrix(0, nrow = length(y), ncol = k+1)
		G.mat[,1]<-0
		G.mat[,k+1]<-1
		for (i in 2:k) { G.mat[,i]<-pnorm(alpha[i-1]+Xb) }
    	pi <- matrix(0, nrow = length(y), ncol = k)
    	for (i in 2:(k + 1)) {
        	pi[, i - 1] <- G.mat[,i] - G.mat[,i-1]
    	}
#		pi <- .Call("do_row_product_sums", pi, Ymat)
#		loglik <- .Call("do_vector_sum", log(pi))
		pi <- pi*Ymat		
		pi <- .Call("do_row_products", pi)
		loglik <- .Call("do_vector_sum", log(pi))		
#		loglik <- .Call("do_vector_sum", log(pi)*Ymat)
	} else if (link=="cloglog") {
		z <- matrix(ncol=k-1, nrow=dim(w)[1])
    	for (i in 1:(k-1)) {
	    	z[,i] <- alpha[i] + Xb
		}
   		pi.z <- matrix(ncol=k, nrow=dim(w)[1])
		for (i in 1:k) {
			if (i==1) {
				pi.z[,i]<-1 - exp(-exp(z[,i]))
			} else if (i <= k-1) {
				pi.z[,i]<-exp(-exp(z[,i-1])) - exp(-exp(z[,i]))
			} else if (i==k) {
				pi.z[,i]<-exp(-exp(z[,i-1]))
			}
		}
#    	pi <- .Call("do_row_product_sums", pi.z, Ymat)
#    	loglik <- .Call("do_vector_sum", log(pi))
		pi <- pi.z*Ymat		
		pi <- .Call("do_row_products", pi)
		loglik <- .Call("do_vector_sum", log(pi))	
#		loglik <- .Call("do_vector_sum", log(pi.z)*Ymat)
	}
    -loglik
}