
##############################################
##   Delsol, Ferraty and Vieu test for no   ##
##      functional-scalar interaction       ##
##############################################

##############################################
## File created by Eduardo Garcia-Portugues ##
## using code from library fda.usc          ##
##############################################


# Delsol, Ferraty and Vieu test statisitc for the simple hypothesis of no interaction
dfv.statistic=function(X.fdata,Y,h=quantile(x=metric.lp(X.fdata),probs=c(0.05,0.10,0.15,0.25,0.50)),K=function(x)2*dnorm(abs(x)),weights=rep(1,dim(X.fdata$data)[1]),d=metric.lp,dist=NULL){
	
	# Check if it is necesary to compute the distances matrix
	if(is.null(dist)) dist=d(X.fdata,X.fdata)
	
	# Statistic
	res=sapply(h,function(hh) mean(drop(K(dist/hh)%*%Y)^2*weights))

	return(res)

}

# Delsol, Ferraty and Vieu test for the simple hypothesis of no interaction with bootstrap calibration
dfv.test=function(X.fdata,Y,B=5000,h=quantile(x=metric.lp(X.fdata),probs=c(0.05,0.10,0.15,0.25,0.50)),K=function(x)2*dnorm(abs(x)),weights=rep(1,dim(X.fdata$data)[1]),d=metric.lp,verbose=TRUE){
	
	# REAL WORLD
	dist=d(X.fdata,X.fdata)
	Tn=dfv.statistic(X.fdata=X.fdata,Y=Y,h=h,K=K,weights=weights,dist=dist)
	
	# BOOTSTRAP WORLD
	Tn.star=matrix(nrow=B,ncol=length(h))
	colnames(Tn.star)=paste("Tn.boot(",sprintf("h=%.3f",h),")",sep="")
	
	if(verbose) pb=txtProgressBar(style=3)
	for(i in 1:B){
		
		# Bootsrtap version of Tn: perturbation with a centred and unit variance noise
		Tn.star[i,]=dfv.statistic(X.fdata=X.fdata,Y=rwild(Y,"golden"),h=h,K=K,weights=weights,dist=dist)
		
		# Progress bar
		if(verbose) setTxtProgressBar(pb,i/B)
		
	}
	
	# P-value
	pvalue=sapply(1:length(h),function(i) sum(Tn.star[,i]>Tn[i]))/B
	names(pvalue)=paste("p-value(",sprintf("h=%.3f",h),")",sep="")

	# Result: class htest
	names(Tn)=paste("Tn(",sprintf("h=%.3f",h),")",sep="")
	result=list(statistic=Tn,boot.statistics=Tn.star,p.value=pvalue,method="Delsol, Ferraty and Vieu test for no functional-scalar interaction",B=B,h=h,K=K,weights=weights,d=d,data.name="Y=0+e")
	class(result)="htest"
	return(result)
	
}
