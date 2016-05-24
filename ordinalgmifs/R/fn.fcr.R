fn.fcr <-
function(par, w, x, beta, y, k, Ymat, Cum.Ymat, link) {
	if (dim(w)[2]!=0) {
    	zeta <- par[k:length(par)]
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
		G.mat <- matrix(0, nrow=length(y), ncol=k-1)  
		for (i in 2:k) {
			G.mat[, i-1] <- G(par[i-1] + Xb, link)
		}
		pi <- Ymat[, 1:(k-1)] * G.mat + (1 - G.mat) * 
			Cum.Ymat * (1 - Ymat[, 1:(k-1)])
		pi <- .Call("do_row_products", pi)
		loglik <- .Call("do_vector_sum", log(pi))
	} else if (link=="probit") {
	    G.mat <- matrix(0, nrow = length(y), ncol = k)  
    	G.mat[,1]<-pnorm(par[1]+Xb)
    	G.mat[,2]<-pnorm(par[2]+Xb)*(1-G.mat[,1])
    	if (k>3) {
    	for (i in 3:(k-1)) {
			Gmatsum <- .Call("do_matrix_sum_rows", G.mat[,1:(i-1)])
			G.mat[,i] <- pnorm(par[i] + Xb) * (1 - matrix(Gmatsum, nrow=nrow(G.mat), byrow=T))
			}
		}
		Gmatsum <- .Call("do_matrix_sum_rows", G.mat[,1:(k-1)])
		G.mat[,k] <- 1 - matrix(Gmatsum, nrow=nrow(G.mat), byrow=T)
    	pi <- Ymat*G.mat
		pi <- .Call("do_row_products", pi)
		loglik <- .Call("do_vector_sum", log(pi))
    } else if (link=="cloglog") {
    	G.mat <- matrix(0, nrow = length(y), ncol = k)  
    	G.mat[,1]<-G(par[1]+Xb, link)
    	G.mat[,2]<-G(par[2]+Xb, link)*(1-G.mat[,1])
    	if (k>3) {
    		for (i in 3:(k-1)) {
    			G.mat[,i]<-G(par[i]+Xb, link)*(1-matrix(apply(G.mat[,1:(i-1)],1,sum),nrow=nrow(G.mat),byrow=T))
			}
		}
		G.mat[,k]<-1-matrix(apply(G.mat[,1:(k-1)],1,sum),nrow=nrow(G.mat),byrow=T)
#		pi <- Ymat*G.mat
#    	pi <- apply(pi,1,sum)
#    	loglik <- sum(log(pi))
		pi <- .Call("do_matrix_sum_rows", Ymat*G.mat)
		loglik <- .Call("do_vector_sum", log(pi))
    }
    -loglik
}
