resid.mix <- function(object,x,y,std=FALSE)
{
	int <- object$intercept
	if(int) x <- cbind(1,x)
	theta <- object$theta
	K     <- length(theta)
	resid <- list()
	gamma <- gfun(x,y,theta)$gamma
	for(k in 1:K) {
		div <- if(std)
			sqrt(hatfun(x,gamma[,k])*theta[[k]]$sigsq)
		else 
			1
		resid[[k]] <- drop(y-x%*%theta[[k]]$beta)/div
	}
	rslt <- list(resid=matrix(unlist(resid),ncol=K),gamma=gamma,
                     x=if(int) drop(x[,-1]) else x,y=y)
	class(rslt) <- "mresid"
	rslt
}
