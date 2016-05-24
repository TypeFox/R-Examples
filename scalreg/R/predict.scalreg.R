predict.scalreg <-
function(object, newX=NULL,...) {
	if(is.null(newX))
		y = fitted(object)
	else{
		y=as.vector(newX%*%object$coefficients)
	}
	y
}
