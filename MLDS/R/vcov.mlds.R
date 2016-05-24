`vcov.mlds`<- function(object, ...) {
	if (object$method == "glm"){
		return(vcov(object$obj))
		} else {
		warning("Fitted values were constrained which may distort standard error estimates!")
		return(solve(object$hess))	
		}
	}
	
`vcov.mlbs`<- function(object, ...) 
		return(vcov(object$obj))
