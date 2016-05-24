CrossValidationMuFicokm <- function(model,indice){
	nlevel <- model$nlevel
	model1 <- model$cok[[nlevel]]
	nugget <- model$nuggets[nlevel]
	aux1 <- covMatrix(model1@covariance, X = model1@X, noise.var = model1@noise.var)
    	C1 <- aux1$C
	nC1 <- dim(C1)[1]
	I1 <- diag(rep(nugget,nC1))
	C1 <- C1 + I1
	iC1 <- solve(C1)

	M1 <- iC1%*%( model1@y -  model1@F%*%model1@trend.coef)
	
	err1 <- solve(iC1[indice,indice],M1[indice,])
	Cov1 <- solve(iC1[indice,indice])
	var1 <- diag(Cov1)

	return(list(CVerr=err1,CVvar=var1,CVCov=Cov1))
}

