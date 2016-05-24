out <- function(x, rho, gamma, scores=FALSE, df.method="reparametrization"){
	if(class(x)!="fanc") stop('the class of object "x" must be "fanc"')
	 if(scores==TRUE && is.null(x$x)==TRUE) stop("Data matrix is needed for computing the factor score in fitting procedure by fanc")
	 
	 gamma_vec <- x$gamma
	 gamma_length <- length(gamma_vec)
	 if(gamma==Inf) gamma_index <- 1
	 if(gamma!=Inf) gamma_index <- which.min(abs(gamma-gamma_vec))
	 rhovec <- x$rho[,gamma_index]
	 rho_index <- which.min(abs(rho-rhovec))	
	 Lambda <- x$loadings[[gamma_index]][[rho_index]]
	 diagPsi <- x$uniquenesses[rho_index,,gamma_index]
	 if(x$cor.factor==TRUE){
		Phi <- x$Phi[,,rho_index,gamma_index]
		Phi <- as.matrix(Phi)
	 }
	 if(scores==TRUE){
		Lambda_mat <- as.matrix(Lambda)
		 diagPsiinvrep <- matrix(diagPsi^(-1),nrow(Lambda),ncol(Lambda))
		 diagPsiinvLambda <- diagPsiinvrep * Lambda_mat
		 M0 <- crossprod(Lambda_mat,diagPsiinvLambda)
		 if(x$cor.factor==TRUE) M <- M0 + solve(Phi)
		 if(x$cor.factor==FALSE) M <- M0 + diag(x$factors)
		 solveM <- solve(M)
		 PsiinvLambdaMinv <-diagPsiinvLambda %*% solveM
		 ans_scores <- x$x %*% PsiinvLambdaMinv
	 }
	if(df.method=="reparametrization") df <- x$df[rho_index,gamma_index]
	if(df.method=="active") df <- x$dfnonzero[rho_index,gamma_index]
	
	
	 if(is.null(x$AIC)==FALSE){
		if(df.method=="reparametrization"){
			AIC <- x$AIC[rho_index,gamma_index]
			BIC <- x$BIC[rho_index,gamma_index]
			CAIC <- x$CAIC[rho_index,gamma_index]
			EBIC <- x$EBIC[rho_index,gamma_index]
		}
		if(df.method=="active"){
			AIC <- x$AIC_dfnonzero[rho_index,gamma_index]
			BIC <- x$BIC_dfnonzero[rho_index,gamma_index]
			CAIC <- x$CAIC_dfnonzero[rho_index,gamma_index]
			EBIC <- x$EBIC_dfnonzero[rho_index,gamma_index]
		}
		 criteria <- c(AIC,BIC,CAIC,EBIC)
		 names(criteria) <- c("AIC","BIC","CAIC","EBIC")
	 }
	 
	 gamma0 <- gamma_vec[gamma_index]
	 rho0 <- rhovec[rho_index]
	 
	 
	if(is.null(x$GFI)==FALSE){
		if(df.method=="reparametrization"){
			GFI <- x$GFI[rho_index,gamma_index];
			AGFI <- x$AGFI[rho_index,gamma_index];
			CFI <- x$CFI[rho_index,gamma_index];
			RMSEA <- x$RMSEA[rho_index,gamma_index];
			SRMR <- x$SRMR[rho_index,gamma_index];
		}
		if(df.method=="active"){
			GFI <- x$GFI[rho_index,gamma_index];
			AGFI <- x$AGFI_dfnonzero[rho_index,gamma_index];
			CFI <- x$CFI_dfnonzero[rho_index,gamma_index];
			RMSEA <- x$RMSEA_dfnonzero[rho_index,gamma_index];
			SRMR <- x$SRMR[rho_index,gamma_index];
		}
		 GOF <- c(GFI,AGFI,CFI,RMSEA,SRMR)
		 names(GOF) <- c("GFI","AGFI","CFI","RMSEA","SRMR")
	 }

	 
	 ans <- list(loadings=Lambda, uniquenesses=diagPsi)
	 if(x$cor.factor==TRUE) ans <- append(ans,list(Phi=Phi))
	 if(scores==TRUE) ans <- append(ans,list(scores=ans_scores)) 
	 ans <- append(ans,list(df=df)) 
	 if(is.null(x$AIC)==FALSE) ans <- append(ans,list(criteria=criteria))
	 if(is.null(x$GFI)==FALSE) ans <- append(ans,list(goodness.of.fit=GOF))
	 ans <- append(ans,list(rho=rho0, gamma=gamma0)) 
	 ans
 }
 
 