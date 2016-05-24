## -----------------------------------------------------------------------------
## Fonction limitFunction
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

#' @import Matrix
limitFunction = function(model,kernel=NULL) {

	if(typeof(model)=="S4") {
		fun = function(X) {
			X = t(as.matrix(X))
			p = predict(model,X,type="UK",checkNames=FALSE)
			return(p)
		}
	}
	else {
		if(is.null(kernel)) {
			kernel_ind = model$kernel
			if(kernel_ind==0) {kernel = function(x,y) {t(x)%*%y}}
			if(kernel_ind==1) {kernel = kernlab::polydot(degree=model$degree, scale=model$gamma, offset=model$coef0)}
			if(kernel_ind==2) {kernel = kernlab::rbfdot(sigma=model$gamma)}
			if(kernel_ind==3) {kernel = kernlab::tanhdot(scale=model$gamma, offset=model$coef0)}
		}
	
		fun = function(X) {
			X = t(as.matrix(X))
		 	K = matrix(kernlab::kernelMatrix(kernel,model$SV,X),dim(model$SV)[1],dim(X)[1])
			res = list(mean = c(t(as.matrix(model$coefs))%*%K - model$rho),
			           sd = 1) #list names are mean and sd to keep consistency with predict output in kriging case
			return(res)
		}
	}
	return(fun)
}