scalreg <-
function(X,y=NULL,lam0=NULL,LSE=FALSE){
	X <- as.matrix(X)
	if(!is.null(y)){
		y <- as.numeric(y)
		est <- slassoEst(X, y, lam0)
		est$fitted.values <- as.vector(X %*% est$coefficients)
		est$residuals <- y - est$fitted.values
		est$type = "regression"
		if(LSE==TRUE){
			lse=lse(X,y,indexset=which(est$coefficients!=0))
			est$lse=lse
		}
	}
	if(is.null(y)){
		est <- slassoInv(X,lam0,LSE)
		est$type = "precision matrix"
	}
	est$call <- match.call()
	class(est) <- "scalreg"
	est
}
