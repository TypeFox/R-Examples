slassoInv <-
function(X,lam0=NULL,LSE=F){
	nX=dim(X)[1]; pX=dim(X)[2]
	hsigma=rep(0,pX);Beta=matrix(-1,pX,pX); res=matrix(0,nX,pX)
	if(LSE==T) {
		hsigma.lse=rep(0,pX);Beta.lse=matrix(-1,pX,pX); res.lse=matrix(0,nX,pX)
	} 
	for(j in 1:pX){
	#	cat(j,"\n")
		scalefac=sqrt(colSums(X[,-j]^2)/nX)
		Xj=t(t(X[,-j])/scalefac)
		objtmp = slassoEst(Xj,X[,j],lam0)
		hsigma[j] = objtmp$hsigma
		Beta[-j,j] = objtmp$coefficients/scalefac
		res[,j]=objtmp$residuals
		
		if(LSE==T) {
			lsetmp=lse(Xj,X[,j],indexset=which(objtmp$coefficients!=0))
			hsigma.lse[j] = lsetmp$hsigma
			Beta.lse[-j,j] = lsetmp$coefficients/scalefac
			res.lse[,j]=lsetmp$residuals
		}
	}
	tTheta=diag(hsigma^(-2)); tTheta=-Beta%*%tTheta
	hTheta=tTheta*(abs(tTheta)<=abs(t(tTheta)))+t(tTheta)*(abs(tTheta)>abs(t(tTheta)))
	est=list(precision=hTheta,hsigma=hsigma)
	
	if(LSE==T) {
		tTheta.lse=diag(hsigma.lse^(-2)); tTheta.lse=-Beta.lse%*%tTheta.lse
		hTheta.lse=tTheta.lse*(abs(tTheta.lse)<=abs(t(tTheta.lse)))+t(tTheta.lse)*(abs(tTheta.lse)>abs(t(tTheta.lse)))
		lse=list(precision=hTheta.lse,hsigma=hsigma.lse)
		est$lse=lse
	}
	
	est
}
