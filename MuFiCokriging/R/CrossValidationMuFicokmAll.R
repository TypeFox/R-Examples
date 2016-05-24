CrossValidationMuFicokmAll <- function(model,indice){

	Nestdesign <- model$Dnest
	nlevel <- model$nlevel
	d <- dim(as.matrix(Nestdesign$PX))[2]
	n <- length(indice)

	indice.nlev <- list()
	for(i in 2:nlevel){
		grid <- expand.grid(as.matrix(ExtractNestDesign(Nestdesign,i-1))[,1],as.matrix(ExtractNestDesign(Nestdesign,nlevel))[indice,1])	
		dist <- (grid[,1]-grid[,2])^2
		if(d > 1){
			for(j in 2:d){
				grid <- expand.grid(as.matrix(ExtractNestDesign(Nestdesign,i-1))[,j],as.matrix(ExtractNestDesign(Nestdesign,nlevel))[indice,j])	
				dist <- dist + (grid[,1]-grid[,2])^2			
			}
		}
		matdist <- matrix(dist,ncol=n)
		indice.nlev[[i-1]] <- max.col(-t(matdist))
	}

	indice1 <- indice.nlev[[1]]
	model1 <- model$cok[[1]]
	aux1 <- covMatrix(model1@covariance, X = model1@X, noise.var = model1@noise.var)
    	C1 <- aux1$C

	nugget1 <- model$nuggets[1]
	nC1 <- dim(C1)[1]
	I1 <- diag(rep(nugget1,nC1))

	iC1 <- solve(C1+I1)

	M1 <- iC1%*%( model1@y -  model1@F%*%model1@trend.coef)
	
	err1 <- solve(iC1[indice1,indice1],M1[indice1,])
	Cov1 <- solve(iC1[indice1,indice1])
	var1 <- diag(Cov1)

	CVerrall <- list()
	CVvarall <- list()
	CVCovall <- list()

	CVerrall[[1]] <- err1
	CVvarall[[1]] <- var1
	CVCovall[[1]] <- Cov1

	indice.nlev[[nlevel]] <- indice

	for(i in 2:nlevel){
		indicei <- indice.nlev[[i]]
		modeli <- model$cok[[i]]
		auxi <- covMatrix(modeli@covariance, X = modeli@X, noise.var = modeli@noise.var)
    		Ci <- auxi$C

		nuggeti <- model$nuggets[i]
		nCi <- dim(Ci)[1]
		Ii <- diag(rep(nuggeti,nCi))

		iCi <- solve(Ci+Ii)

		Mi <- iCi%*%( modeli@y -  modeli@F%*%modeli@trend.coef)
	
		erri <- solve(iCi[indicei,indicei],Mi[indicei,])
		Covi <- solve(iCi[indicei,indicei])
		vari <- diag(Covi)

		rhoi <- as.numeric(modeli@trend.coef[1])

		erri <- solve(iCi[indicei,indicei],Mi[indicei,])
		Covi <- solve(iCi[indicei,indicei])
		vari <- diag(Covi)

		CVerrall[[i]] <- rhoi*CVerrall[[i-1]] + erri
		CVvarall[[i]] <- rhoi^2*CVvarall[[i-1]] + vari
		CVCovall[[i]] <- rhoi^2*CVCovall[[i-1]] + Covi
	
	}
	

	return(list(CVerrall=CVerrall[[nlevel]],CVvarall=CVvarall[[nlevel]],CVCovall=CVCovall[[nlevel]], CVerr=CVerrall, CVvar=CVvarall, CVCov=CVCovall))
}
