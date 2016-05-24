FUN.iterate.GCV <- function(TR.Y.mat, TR.X.mat, N, T, P){
	y     <- matrix(TR.Y.mat,  nrow= N*T, ncol = 1)			# (TN x 1)
	x     <- matrix(TR.X.mat,  nrow= N*T, ncol = P)			# (TN x P)
	t.seq <- seq(0, 1, length.out=T)

	xx 	   <- crossprod(x)				# (p x p)
	inv.xx   <- solve(xx)					# (p x p)
	inv.xx.x <- inv.xx%*%t(x)				# (p x TN)

  ## function that will be iterated 
	FUN.ols.beta <- function(updated.y, x, inv.xx.x){
		beta <- tcrossprod(inv.xx.x, t(updated.y))
		}


  ## Starting values 		
	ymats <- matrix(y, T, N)
	xmats <- matrix(x, T, (N*P))
	zmats <- cbind(ymats, xmats)
	trzma <- zmats - svd.pca(zmats, given.d = round(sqrt(min(N, T))) )$Q.fit
	tryxm <- matrix(trzma, (T*N), (P+1))
	try   <- tryxm[, 1, drop = FALSE]
	trx   <- tryxm[, -1, drop = FALSE]
	beta.0 <- coef(lm(try ~ -1 + trx ))			# (p x p)  


  ## Iteration
	inner.iteration <- function(y, x, inv.xx.x =inv.xx.x, 
				  beta.0 = beta.0, t.seq,i=1){
	nr <- length(t.seq)
	nc <- nrow(y)/nr
  ## Iteration (0): initial cumputations
	## w.0 = y-x%*%beta.0 
		w.0 <- y - tcrossprod(x, t(beta.0))

  	## W.0: write w.0 in a matrix form 
		W.0 <- matrix(w.0, nr, nc)

 	## PCA.0: PCA.0 computation, OptDim.0 and y.fitted
		PCA.0    <- smooth.Pspline(x = t.seq, y = W.0, method = 3)

		y.fitted.0 = PCA.0$ysmth

   ## Iteration (+1)
  	## y.updated.0: updat y.updated.0 = y - fs.0
		y.updated.1 <- y -  c(y.fitted.0)

  	## beta.1: OLS.1 computation for the computed fs.0 in interation 0 
		beta.1 <- FUN.ols.beta(y.updated.1, x, inv.xx.x) 

  	## convergence condition
		if(all( abs((beta.0 - beta.1)) < 1e-3)| i  == 100){
			Result <- list(PCA=PCA.0, beta=beta.1, Nbr.Iterations = i)
			Result
			}

		else inner.iteration(y=y, x=x, inv.xx.x =inv.xx.x 
			,beta.0 = beta.1, t.seq = t.seq, (i+1))
			
	}
	
	Result <- inner.iteration(y = y, x= x, t.seq = t.seq, 
					inv.xx.x = inv.xx.x, 
					beta.0 = beta.0, i=1)

}

