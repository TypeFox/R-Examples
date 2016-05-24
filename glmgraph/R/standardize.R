standardizeX <- function(X){
	n <- nrow(X)
	p <- ncol(X)
	Xc <- apply(X,2,mean)
	X  <- sweep(X,2,Xc,FUN="-")
	Xs <- apply(X,2,sd)
	X  <- sweep(X,2,Xs,FUN="/")
	return(list(X,Xc,Xs))
}








