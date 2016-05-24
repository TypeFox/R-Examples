## ' Similar to \code{\link[KrigInv]{computeQuickKrigcov}}, but without dataframe to matrix conversion
## ' for faster execution
## ' @title Quick computation of kriging covariances
## ' @param model A Kriging model of \code{\link[DiceKriging]{km}} class.
## '  @param integration.points p*d matrix of points for numerical integration in the X space.
## ' @param X.new The new point where we calculate kriging covariances. The calculated covariances 
## ' are the covariances between this new point and all the integration points.
## ' @param precalc.data List containing precalculated data. This list is generated using the function 
## ' \code{\link[KrigInv]{precomputeUpdateData}}
## ' @param F.newdata The value of the kriging trend basis function at point X.new
## ' @param c.newdata The (unconditional) covariance between X.new and the design points
## ' @details See \code{\link[KrigInv]{computeQuickKrigcov}} for more details
## ' 
## ' @return \code{TRUE} if the point should not be tested
## ' @export

computeQuickKrigcov2 <- function(model,integration.points,X.new, precalc.data, F.newdata , c.newdata){
#   X.new <- as.matrix(X.new)
#   if (model@d == 1){ integration.points <- as.matrix(integration.points) } 
#   colnames(integration.points) <- colnames(model@X)
	c.xnew.integpoints <- covMat1Mat2(X1=integration.points,X2=X.new, object=model@covariance, nugget.flag=model@covariance@nugget.flag)
	cov.std <- c.xnew.integpoints - crossprod(precalc.data$Kinv.c.olddata,c.newdata)		#O(M.n) here
  if (is.null(F.newdata))
  {kn=cov.std
  } else
  { second.member <- t(F.newdata - crossprod(c.newdata,precalc.data$Kinv.F))
    cov.F <- precalc.data$first.member%*%second.member
    kn <- cov.F+cov.std}
	return(kn)
}

