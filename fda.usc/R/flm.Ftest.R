
##############################################
##    Functional F-test (beta=0) in FLM     ##
##############################################

##############################################
## File created by Eduardo Garcia-Portugues ##
## using code from library fda.usc          ##
##############################################


# FLM F-Test statistic
Ftest.statistic=function(X.fdata,Y){
	
	# Statistic
	res=as.numeric(norm.fdata(func.mean(fdata.cen(X.fdata)$Xcen*(Y-mean(Y)))))

	return(res)
	
}

# FLM F-Test with bootstrap calibration
flm.Ftest=function(X.fdata,Y,B=5000,verbose=TRUE){
	
	# REAL WORLD
	Tn=Ftest.statistic(X.fdata=X.fdata,Y=Y)
	
	# BOOTSTRAP WORLD
	Tn.star=numeric(B)
	if(verbose) pb=txtProgressBar(style=3)
	for(i in 1:B){
		
		# Bootsrtap version of Tn: perturbation with a centred and unit variance noise
		Tn.star[i]=Ftest.statistic(X.fdata=X.fdata,Y=rwild(Y,"golden"))
		
		# Progress bar
		if(verbose) setTxtProgressBar(pb,i/B)
		
	}
	
	# P-value
	pvalue=sum(Tn.star>Tn)/B
	
	# Result: class htest
	names(Tn)="F-test"
	result=structure(list(statistic=Tn,boot.statistics=Tn.star,p.value=pvalue,method="Functional Linear Model F-test",B=B,data.name="Y=<X,0>+e"))
	class(result)="htest"
	return(result)
	
}
