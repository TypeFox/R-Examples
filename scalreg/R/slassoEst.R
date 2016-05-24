slassoEst <-
function(X,y,lam0=NULL){
	nX=dim(X)[1]; pX=dim(X)[2]
	if(is.null(lam0)){
		if(pX>10^6){
			lam0 = "univ"
		} else lam0 = "quantile"
	} 
	if(lam0=="univ" | lam0=="universal")
		lam0=sqrt(2*log(pX)/nX)
	if(lam0=="quantile"){
		L=0.1; Lold=0
		while(abs(L-Lold)>0.001){
			k=(L^4+2*L^2); Lold=L; L=-qnorm(min(k/pX,0.99)); L=(L+Lold)/2
		}
		if(pX==1) L=0.5
		lam0 = sqrt(2/nX)*L
	}

	objlasso=lars(X,y,type="lasso",intercept=FALSE, normalize=FALSE, use.Gram=FALSE)
#	cat("Lasso path created!\ttime = ",date(),"\n")
	
	sigmaint=0.1; sigmanew=5; flag=0
	while(abs(sigmaint-sigmanew)>0.0001 & flag <= 100){
		flag=flag+1
		sigmaint=sigmanew; lam=lam0*sigmaint
		hy=predict.lars(objlasso,X,s=lam*nX,type="fit", mode="lambda")$fit
		sigmanew=sqrt(mean((y-hy)^2))
	}
	hsigma=sigmanew; hlam=lam
	hbeta=predict.lars(objlasso,X,s=lam*nX,type="coefficients", mode="lambda")$coef
	hy=predict.lars(objlasso,X,s=lam*nX,type="fit", mode="lambda")$fit
	return(list(hsigma=hsigma,coefficients=hbeta,residuals=y-hy, fitted.values=hy))
}
