glmnetEmpty <-
function(X, lambda){
	fit = list()
	
	fit$a0 <- rep(0, length(lambda))
	fit$lambda <- lambda
	fit$df <- 0
	fit$dim <- c(ncol(X), length(lambda))
	
	fit$beta <- Matrix::Matrix(0, ncol(X), length(lambda))
	rownames(fit$beta) <- colnames(X)
	colnames(fit$beta) <- paste("s", 0:(length(lambda)-1), sep="")
	
	return(fit)
}
