precomputeUpdateData <- function(model,integration.points){
	
	#precalculates the covariances between the integration points and the design points
  integration.points <- as.matrix(integration.points)
  colnames(integration.points) <- colnames(model@X)
  
	c.olddata <- covMat1Mat2(X1=model@X, X2=integration.points, object=model@covariance, nugget.flag=model@covariance@nugget.flag)	
	Tinv.c.olddata <- backsolve(t(model@T), c.olddata, upper.tri=FALSE)						#a big backsolve
	Kinv.c.olddata <- backsolve(model@T, Tinv.c.olddata, upper.tri=TRUE)					#this
	Kinv.F <- backsolve(model@T,model@M)													#this
	
	inv.tF.Kinv.F<- chol2inv(chol(crossprod(model@M)))
	f.integration.points <- model.matrix(model@trend.formula, data=data.frame(integration.points))
	first.member <- (f.integration.points-crossprod(c.olddata,Kinv.F))%*%inv.tF.Kinv.F		#this
	
	return(list(
	Kinv.c.olddata=Kinv.c.olddata,
	Kinv.F=Kinv.F,
	first.member=first.member
	))
}

