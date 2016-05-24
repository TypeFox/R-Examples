bothsidesmodel.lrt <-
function(x,y,z,pattern0,patternA=matrix(1,nrow=ncol(x),ncol=ncol(z))) {
	bsmA <- bothsidesmodel.mle(x,y,z,patternA)
	bsm0 <- bothsidesmodel.mle(x,y,z,pattern0)
    chisq <- bsm0$Deviance-bsmA$Deviance
    df <- bsmA$Dim - bsm0$Dim
	list(Chisq=chisq,df=df,pvalue=1-pchisq(chisq,df))
}
