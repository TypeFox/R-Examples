dlogr.step <-
function (x.t,y.t,betahat.tm1,varbetahat.tm1,tune.mat) {

	if (!is.matrix(x.t)) {
	   dim(x.t) <- c(1,length(x.t))
	}
    # first decide the forgetting factor
    temp <- apply(tune.mat,1,laplace.fn,x.t=x.t,y.t=y.t,betahat.tm1=betahat.tm1,varbetahat.tm1=varbetahat.tm1)
    lambda <- tune.mat[which.max(temp),]
    Rhat.t <- varbetahat.tm1
    diag(Rhat.t) <- diag(Rhat.t) / lambda
    laplace.t=max(temp)
    yhat.t <- dlogr.predict(x.t,betahat.tm1)
   
    Del1 <- t(x.t) %*% (y.t - yhat.t)
	Del2 <- -solve(Rhat.t) - (t(x.t) * matrix(rep(yhat.t*(1-yhat.t),dim(x.t)[2]),nrow=dim(x.t)[2],byrow=TRUE)) %*% x.t

    betahat.t <- betahat.tm1 - (solve(Del2) %*% Del1)
    varbetahat.t <- solve(-Del2)
	  diag(varbetahat.t) <- abs(diag(varbetahat.t))
  	return(list(betahat.t=betahat.t,varbetahat.t=varbetahat.t,laplace.t=laplace.t))
}

