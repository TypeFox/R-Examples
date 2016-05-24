confintAsin.fun<-function(mlePP, level=0.95)
{
	coef<-mlePP@coef
	if (dim(mlePP@vcov)[1]==0)  stop('The confidence intervals cannot be calcualted since the
covariance matrix of the coefficientes cannot be estimated')
	sterr<-diag(mlePP@vcov)**0.5
	CIupper<-coef+sterr*qnorm(1-(1-level)/2)
	CIlower<-coef-sterr*qnorm(1-(1-level)/2)
	CI<-cbind(CIlower, CIupper)
	a <- (1 - level)/2
 	a <- c(a, 1 - a)
    	pct <- paste(round(100 * a, 1), "%")
	dimnames(CI)<-list(names(coef), pct)
	return(CI)
}

