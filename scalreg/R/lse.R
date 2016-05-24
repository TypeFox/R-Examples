lse <-
function(X,y,indexset){
	hbeta=rep(0,dim(X)[2])
	if(length(indexset)>0){
		objlm=lm(y~X[,indexset]+0)
		hbeta[indexset]=objlm$coefficients
		hsigma=sqrt(mean(objlm$residuals^2))
		residuals=objlm$residuals; fitted.values=objlm$fitted.values
	} else{
		hsigma=sqrt(mean(y^2)); residuals=y; fitted.values=rep(0,dim(X)[2])
	}
	return(list(hsigma=hsigma,coefficients=hbeta,residuals=residuals, fitted.values=fitted.values))
}
