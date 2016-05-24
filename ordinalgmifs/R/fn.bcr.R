fn.bcr <-
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
    	G.mat <- matrix(0, nrow=nrow(Ymat), ncol = k)  
		G.mat[,k] <- G(par[k-1] + Xb, link)
		G.mat[,k-1] <- G(par[k-2]+Xb, link) * (1 - G.mat[,k])
		if (k > 3) {
			for (i in (k-2):2) {
				Gmatsum <- .Call("do_matrix_sum_rows", G.mat[,k:(i+1)])
				G.mat[,i] <- G(par[i-1] + Xb, link) * (1 - matrix(Gmatsum, nrow=nrow(G.mat), byrow=T))
			}
		}
		Gmatsum <- .Call("do_matrix_sum_rows", G.mat[,k:2])
		G.mat[,1] <- 1 - matrix(Gmatsum, nrow=nrow(G.mat), byrow=T)
		pi <- Ymat * G.mat
		pi <- .Call("do_matrix_sum_rows", pi)
		loglik <- .Call("do_vector_sum", log(pi))
	} else if (link=="probit") {
	    G.mat <- matrix(0, nrow = length(y), ncol = k)  
    	G.mat[,k]<-pnorm(par[k-1]+Xb)
    	G.mat[,k-1]<-pnorm(par[k-2]+Xb)*(1-G.mat[,k])
    	if (k>3) {
    		for (i in (k-2):2) {
				Gmatsum <- .Call("do_matrix_sum_rows", G.mat[,k:(i+1)])
				G.mat[,i] <- pnorm(par[i-1] + Xb) * (1 - matrix(Gmatsum, nrow=nrow(G.mat), byrow=T))
			}
		}
		Gmatsum <- .Call("do_matrix_sum_rows", G.mat[,k:2])
		G.mat[,1] <- 1 - matrix(Gmatsum, nrow=nrow(G.mat), byrow=T)
    	pi <- Ymat*G.mat
 		pi <- .Call("do_matrix_sum_rows", pi)
		loglik <- .Call("do_vector_sum", log(pi))
	} else if (link=="cloglog") {
    	G.mat <- matrix(0, nrow = length(y), ncol = k)  
    	G.mat[,k]<-G(par[k-1]+Xb, link)
    	G.mat[,k-1]<-G(par[k-2]+Xb, link)*(1-G.mat[,k])
    	if (k>3) {
    		for (i in (k-2):2) {
    			G.mat[,i]<-G(par[i-1]+Xb, link)*(1-matrix(apply(G.mat[,k:(i+1)],1,sum),nrow=nrow(G.mat),byrow=T))
			}
		}
		G.mat[,1]<-1-matrix(apply(G.mat[,k:2],1,sum),nrow=nrow(G.mat),byrow=T)
#    	pi <- Ymat*G.mat
#    	pi <- rowSums(pi,na.rm=TRUE)
#    	loglik <- .Call("do_vector_sum", log(pi))
		pi <- .Call("do_matrix_sum_rows", Ymat*G.mat)
		loglik <- .Call("do_vector_sum", log(pi))
	}
	-loglik
}