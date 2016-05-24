coef.qut <- 
function(object, mode='glm',...){
	if(mode=='lasso') return(object$beta)
	else{
		if(!attr(object$betaglm,'converged')){
			warning('GLM did not converged, returning LASSO coefficients')
			return(object$lasso)
		}
		else return(object$betaglm)
	}
}