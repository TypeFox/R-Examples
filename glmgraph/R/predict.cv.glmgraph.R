predict.cv.glmgraph <- function(object,X,s=c("lambda1.min","lambda1.1se"),type=c("response", "coefficients", "class", "nzeros","link"),...) {
	s <- match.arg(s)
	if(s=="lambda1.min") return(predict(object$obj.min,X,type=type)[[1]])
	else if(s=="lambda1.1se") return(predict(object$obj.1se,X,type=type)[[1]])
	else stop("Invalid type")
}
