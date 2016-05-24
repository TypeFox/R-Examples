# last modified 2011-08-06 by J. Fox

fscores <- function(model, ...){
	UseMethod("fscores")
}

fscores.sem <- function(model, data=model$data, center=TRUE, scale=FALSE, ...){
	m <- model$m
	P <- model$P
	A <- model$A
	var.names <- model$var.names
	observed <- var.names %in% rownames(model$C)
	if (all(observed)) stop("there are no latent variables")
	IAinv <- solve(diag(m) - A)
	Sigma <- IAinv %*% P %*% t(IAinv)
	B <- solve(Sigma[observed, observed]) %*% Sigma[observed, !observed]
	rownames(B) <- var.names[observed]
	colnames(B) <- var.names[!observed]
	if (is.null(data)) return(B)
	X <- as.matrix(data[,var.names[observed]])
	if (center || scale) X <- scale(X, center=center, scale=scale)
	X %*% B
}