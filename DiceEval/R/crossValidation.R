crossValidation <- function(model,K){
	X <- model$data$X	# design of experiments
	Y <- model$data$Y	# input values at the points of the design X
	p <- dim(X)[2]
	n <- dim(X)[1]
	type <- model$type
	model_ <- model$model 	# Model fitted with all observations = model
	classe <- foldsComposition(n,K)
	
	if (K < 2) {
        stop("The number of folds (argument \'K\') should be greater than or equal to 2")
    } else if (K > n) {
        stop("K should be less than or equal to the number of observations")
    }

	
	if(type=="Kriging"){
		theta <- matrix(0,nrow=K,ncol=p)
		trend <- matrix(0,nrow=K,ncol=length(model_@trend.coef))
		if(model_@covariance@param.n==2*p){
			shape <- matrix(0,nrow=K,ncol=p)
		}
	}
	
	Ypred  <- vector("numeric",length=n)
	RMSE_CV <- vector("numeric",length=K)
	MAE_CV <- vector("numeric",length=K)
	for (i in 1:K){
		indice 		<- classe==i
		design 		<- X[!indice,]
		Ycv			<- Y[!indice]
		if (type=="Kriging"){
			# fitting the 'km' model with the same parameter as the input 'km' model		
			modi <- modelFit(X[!indice,],Y[!indice],type=type,formula=model_@trend.formula,
				covtype=model_@covariance@name, nugget = model_@covariance@nugget, 
				nugget.estim = model_@covariance@nugget.estim, noise.var = model_@noise.var, 
				penalty = model_@penalty, optim.method = model_@optim.method, 
				lower = model_@lower, upper = model_@upper, parinit = model_@parinit,
				control = model_@control, gr = model_@gr)
				# ---- tous les modèles sont estimés avec la même initialisation (parinit ?)
				theta[i,] <- modi$model@covariance@range.val
				trend[i,] <- modi$model@trend.coef	
				if(model_@covariance@param.n==2*p){
					shape[i,] <- modi$model@covariance@shape.val
				}
		} else {
			fmla=NA; degree <- NA; a <- NA; k <- NA	
			switch(type,
				"Linear"     = { fmla <- model$formula},
				"Additive"   = { fmla <- model$formula},
				"StepLinear" = { fmla <- model$formula; k <- model$penalty},
				"PolyMARS"   = { a <- model$gcv},
				"MARS"       = { degree <- model$degree})
				modi <- modelFit(design,Ycv,type=type,formula=fmla,degree=degree,gcv=a,penalty=k)
		}
		
		# predicting at testdata points
		Ypred[indice] <- modelPredict(modi,newdata=X)[indice]
		RMSE_CV[i] <- RMSE(Y[indice],Ypred[indice])
		MAE_CV[i] <- MAE(Y[indice],Ypred[indice])
	}
	
	folds <- vector("list", K)
	for(k in 1:K){
		folds[[k]] <- seq(1:n)[classe==k]
	}
	
	Q2 = R2(Y,Ypred)
	if (type == "Kriging"){
		if(model_@covariance@param.n==2*p){
			out <- list(Ypred=Ypred,Q2 = Q2,theta=theta,shape=shape,trend=trend,folds=folds,RMSE_CV=mean(RMSE_CV),MAE_CV=mean(MAE_CV))
		} else out <- list(Ypred=Ypred,Q2 = Q2,theta=theta,trend=trend,folds=folds,RMSE_CV=mean(RMSE_CV),MAE_CV=mean(MAE_CV))
	} else out <- list(Ypred=Ypred,Q2=Q2,folds=folds,RMSE_CV=mean(RMSE_CV),MAE_CV=mean(MAE_CV))
	
	
	return(out)
}

#--------------------#
#  FoldsComposition  #
#--------------------#

foldsComposition <- function(n,K){

	if (K==n){
		classe <- seq(1,n,by=1)
	} else {
	  nb_elem	<- rep(floor(n/K),K)
	  a_ajouter	<- n - floor(n/K)*K
	  if (a_ajouter>0){
	    nb_elem[1:a_ajouter]	<- nb_elem[1:a_ajouter]+1
	  }
	  if (sum(nb_elem)!=n) stop("Problem during the construction of the folds.")

	  classe <- vector("numeric",length=n)
	  i <- 1
	  while (i<=n){
		tmp <- floor(runif(1, min=1, max=K+1))
	    if (nb_elem[tmp]!=0){
		classe[i]	 <- tmp
		nb_elem[tmp] <- nb_elem[tmp]-1
		i 		 <- i+1
	    }
	  }
	}
	return(classe)
  }