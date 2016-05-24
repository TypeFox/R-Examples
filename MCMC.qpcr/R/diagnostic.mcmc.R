diagnostic.mcmc <-
function(model,...) {
	par(mfrow=c(2,2))
	pred=predict(model,type="terms")
	lat=as.vector(apply(model$Liab,2,mean))
	res=lat-pred
	ress=res/apply(res,2,sd)
	plot(res~pred,xlab="predicted value",ylab="lognormal residuals",mgp=c(2.3,1,0),main="residuals vs predicted",...)
	lines(lowess(pred,res),col="red")
	plot(sqrt(abs(res))~pred,xlab="predicted value",ylab="sqrt(standardized residuals)",mgp=c(2.3,1,0),main="Scale-Location",...)
	lines(lowess(pred,sqrt(abs(res))),col="red")
	qqnorm(ress,mgp=c(2.3,1,0),ylab="standardized residuals",main="normal QQ plot",...)
	abline(0,1,col="red")
}
