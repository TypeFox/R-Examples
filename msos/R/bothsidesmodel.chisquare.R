bothsidesmodel.chisquare <-
function(x,y,z,pattern0,patternA=matrix(1,nrow=ncol(x),ncol=ncol(z))) {
	bsm <- bothsidesmodel(x,y,z,patternA)
	which <- patternA*(1-pattern0)
	which <- c(t(which)) == 1
	theta <- c(t(bsm$Beta))[which]
	covtheta <- bsm$Covbeta[which,which]
	chisq <- theta%*%solve(covtheta,theta)
	df <- sum(which)
	list(Theta=theta,Covtheta = covtheta,df = df, Chisq=chisq,pvalue=1-pchisq(chisq,df))
}
