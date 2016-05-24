gdf <- function(object){
	if(object$family != "binomial") stop("gdf function is defined only for binomial family")
	y <- object$y
	X <- object$X
	beta <- object$beta
	beta_dim <- dim(beta)
	p <- beta_dim[1]-1
	np <- beta_dim[2]
	gdf_v <- vector(length=np)
	out.glm <- glm(y~X,family="binomial")
	if(!out.glm$converged){
		warning("complexity was set equal to 'df'")
		gdf_v <- colSums(abs(beta)>0)
	} else {
		mu <- out.glm$fit
		V <- binomial()$variance(mu)
		for(i in 1:np){
			bh <- beta[,i,drop=TRUE]
			A <- which(abs(bh[-1])>0)
			if(length(A)==0){
				Xa <- rep(1,length(y))
				eta <- rep(bh[1],length(y))
			} else {
				Xa <- cbind(1,X[,A])
				ba <- bh[c(1,A+1)]
				eta <- drop(tcrossprod(ba,Xa))
			}
			mu_g <- binomial()$linkinv(eta)
			V_g <- binomial()$variance(mu_g)
			Ib <- crossprod(sqrt(V)*Xa)
			invIb_g <- solve(crossprod(sqrt(V_g)*Xa))
			gdf_v[i] <- drop(crossprod(as.vector(invIb_g),as.vector(Ib)))			
		}
	}
	gdf_v
}