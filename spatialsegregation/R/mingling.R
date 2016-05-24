# mingling.R
# 
# Spatial mingling index for multitype spatial point pattern and graphs
# 
# Author: Tuomas Rajala <tarajala@maths.jyu.fi>
###############################################################################

# minglingF


minglingF<-function(X, r=NULL, target=NULL, ratio=FALSE, ...)
#Mingling index for geometric graph with various range-parameters
#If target type is not given give mean over all types 
#If ratio=TRUE  use M = (1 - M)/(lambda_t/lambda) instead.
{
	# check that X is ppp-object
	verifyclass(X, "ppp")
	if(length(levels(X$marks))<2) stop("Use only on a multitype point pattern (data.frame-marks not yet supported).")
	# if no target given, calculate for all types
	if(is.null(target))
	{
		targeti <- 0
		valu   <- "MinglingMean"
	}
	# else convert to an integer
	else
	{
		if(!is.factor(X$marks))warning("Marks of X are not in factor form. Transforming.")
		X$marks<-as.factor(X$marks)
		targeti<- which( levels(X$marks)  == target)
#		targeti<-which( union(X$marks, NULL) == target)
		if(length(targeti)!=1) stop("Target type not one of pattern types.")
	}
	
	# ratio=logical if an intensity weighted version should be calculated
	if(ratio) funtype <- "Intensity weighted-Mingling Index"
	else funtype <- "Mingling index"
	
	# use the main calc function
	res<-segregationFun(X=X, fun="mingling", r, funpars=c(targeti, as.numeric(ratio)), ...)
	
	# theoretical values in CSR
	sum0<-summary(X)
	l<-sum0$int
	if(!ratio)
		if(targeti==0) theo<-mean(1-sum0$marks[,3]/sum0$int)
		else theo<-1-sum0$marks[targeti,3]/sum0$int
	else theo<-1
	
#	omarks<-order(union(X$marks,NULL)) # right order of marks
#	f<-freqs(X)[omarks]
#	w<-f/sum(f)
#	
	# create the fv-object
	mingling.final<-fv(data.frame(theo=theo,par=res$parvec), 
			argu="par",
			alim=range(res$parvec),
			ylab=substitute(MinglingIndex, NULL),
			desc=c("CSR values","Parameter values"),
			valu="theo",
			fmla=".~par",
			unitname=res$unitname,
			fname=funtype
	)
#	return(res)
	# add all typewise values if no target type given
	if(targeti==0)
	{
		# the values from calculation
		tw<-res$v
		
		# set the names right, and don't forget to check inclusion (might drop some types off)
		colnames(tw)<-union(X$marks[res$included],NULL)
		
		mingling.final<-bind.fv(x=mingling.final,
				y=tw,
				desc=paste("Typewise",funtype,"for type",colnames(tw)),
				labl=colnames(tw)
		)
		
		mingling.final<-bind.fv(x=mingling.final,
				y=data.frame("MinglingMean"=apply(res$v,1,mean,na.rm=TRUE)),
				desc=paste("Mean",funtype,"over types"),
				labl="MinglingMean",
				preferred="MinglingMean"
		)		
		# a frequency weighted mean instead of just a mean, w=freqs/sum(freqs)
		#Iw=apply(res$v,2,weighted.mean,w=w,na.rm=TRUE), 
	}
	
	# if target type given add the values for the target type
	else
	{
		mingling.final<-bind.fv(x=mingling.final,
				y=data.frame("Mingling"=res$v[,1]),
				desc=paste(funtype,"for type",target),
				labl="Mingling",
				preferred="Mingling"
		)
	}
	
	# attach the frequencies too
	attr(mingling.final,"frequencies")<-freqs(X[res$included])
	
	# and some notes
	attr(mingling.final,"neighbourhoodType")<-res$ntype
	attr(mingling.final,"note")<-res$note
	
	# point values
    attr(mingling.final,"point.values")<-res$point.values
	
	# return 
	mingling.final
}

###############################################################################


###############################################################################
mingling.index<-function(X, r=4, ntype="knn", ...)
{
	if(length(r)>1)stop("Use minglingF for vector of parameter values.")
	I0<-minglingF(X, r=r, ntype=ntype, ...)
	data.frame(MinglingMean=I0$MinglingMean, par=I0$par)
}
#eof
