select <- function(x, criterion=c("BIC","AIC","CAIC","EBIC"), gamma,  scores=FALSE, df.method="active"){
	if(class(x)!="fanc") stop('the class of object "x" must be "fanc"')
	if(!missing(gamma)){
		if(gamma<=1) stop("gamma must be greater than 1")
	}
	if(scores==TRUE && is.null(x$x)==TRUE) stop("Data matrix is needed for computing the factor score in fitting procedure by fanc")
	if(is.null(x$AIC)==TRUE) stop("The model selection criterion was not able to be calculated. Data matrix or the number of observations is needed in fitting procedure by fanc.")
    cand <- c("BIC", "AIC", "CAIC", "EBIC")
	criterion <- criterion[1]
	if(sum(criterion==cand) != 1) stop('"criterion" must be "AIC", "BIC, "CAIC" or "EBIC".')

	
	if(df.method=="reparametrization"){
		if(criterion=="AIC") criterion_vec <- x$AIC
		if(criterion=="BIC") criterion_vec <- x$BIC
		if(criterion=="CAIC") criterion_vec <- x$CAIC
		if(criterion=="EBIC") criterion_vec <- x$EBIC
	}
	if(df.method=="active"){
		if(criterion=="AIC") criterion_vec <- x$AIC_dfnonzero
		if(criterion=="BIC") criterion_vec <- x$BIC_dfnonzero
		if(criterion=="CAIC") criterion_vec <- x$CAIC_dfnonzero
		if(criterion=="EBIC") criterion_vec <- x$EBIC_dfnonzero
	}
 

	gamma_vec <- x$gamma
	gamma_length <- length(gamma_vec)
	if(missing(gamma)) gamma_index <- which.min(apply(criterion_vec,2,min))	
	else if(gamma==Inf) gamma_index <- 1
	else if(gamma!=Inf) gamma_index <- which.min(abs(gamma-gamma_vec))


	if(gamma_length == 1) criterion_vec2=c(criterion_vec)
	else criterion_vec2=criterion_vec[,gamma_index]
	
	rho_index <- which.min(criterion_vec2)
	Lambda <- x$loadings[[gamma_index]][[rho_index]]
	diagPsi <- x$uniquenesses[rho_index,,gamma_index]
	 if(x$cor.factor==TRUE){
		Phi <- x$Phi[,,rho_index,gamma_index]
		Phi <- as.matrix(Phi)
	 }
	rho0 <- x$rho[rho_index,gamma_index]
	gamma0 <- gamma_vec[gamma_index]
	criterion_minimum <- min(criterion_vec2)
	if(df.method=="reparametrization") df <- x$df[rho_index,gamma_index]
	if(df.method=="active") df <- x$dfnonzero[rho_index,gamma_index]
	
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
	 if(criterion=="AIC") ans <- append(ans,list(AIC=criterion_minimum))
	 if(criterion=="BIC") ans <- append(ans,list(BIC=criterion_minimum))
	 if(criterion=="CAIC") ans <- append(ans,list(CAIC=criterion_minimum))
	 if(criterion=="EBIC") ans <- append(ans,list(EBIC=criterion_minimum))
	 if(is.null(x$GFI)==FALSE) ans <- append(ans,list(goodness.of.fit=GOF))
	 ans <- append(ans,list(rho=rho0, gamma=gamma0)) 
	 

	 ans
 }
 
 