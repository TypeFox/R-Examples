
computeQuickKrigcov <- function(model,integration.points,X.new,
		precalc.data, F.newdata , c.newdata){
	
  integration.points <- as.matrix(integration.points)
  colnames(integration.points) <- colnames(model@X)
  
	c.xnew.integpoints <- covMat1Mat2(X1=integration.points,X2=X.new, object=model@covariance, nugget.flag=model@covariance@nugget.flag)
	second.member <- t(F.newdata - crossprod(c.newdata,precalc.data$Kinv.F))
	cov.F <- precalc.data$first.member%*%second.member
	cov.std <- c.xnew.integpoints - crossprod(precalc.data$Kinv.c.olddata,c.newdata)		#O(M.n) here
	kn <- cov.F+cov.std
	return(kn)
}

