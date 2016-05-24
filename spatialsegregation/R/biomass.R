# biomass.R
# 
# Neighbourhood mass-element sum.
# 
# Author: Tuomas Rajala <tarajala@jyu.fi>
###############################################################################

# biomassF


biomassF<-function(X, r=NULL, target=NULL, v2=FALSE, ...)
{
	# check that X is ppp-object
	verifyclass(X, "ppp")
	if(length(levels(X$marks))<2) warning("Expected  multitype point pattern (data.frame-marks not yet supported).")
	if(length(X$mass)<X$n) stop("Put the biomass information into $mass-element (vector of length n).")
	# if no target given, calculate for all types
	if(is.null(target))
	{
		targeti <- 0
		valu   <- "Biomass of all species"
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
	
	funtype <- "Biomass sum"
	if(v2)funtype<-"Average biomass"
	# use the main calc function
	res<-segregationFun(X=X, fun="biomass", r, funpars=c(targeti,as.integer(v2)), ...)
	
	theo<-ifelse(v2, mean(X$mass),0)

	# create the fv-object
	biomass.final<-fv(data.frame(theo=theo,par=res$parvec), 
			argu="par",
			alim=range(res$parvec),
			ylab=substitute(Biomass, NULL),
			desc=c("CSR values","Parameter values"),
			valu="theo",
			fmla=".~par",
			unitname=res$unitname,
			fname=funtype
	)

	# add all typewise values if no target type given
	if(targeti==0)
	{
		# the values from calculation
		tw<-res$v
		
		# set the names right, and don't forget to check inclusion (might drop some types off)
		colnames(tw)<-union(X$marks[res$included],NULL)
		
		biomass.final<-bind.fv(x=biomass.final,
				y=tw,
				desc=paste("Typewise neighbourhood",funtype,"for type",colnames(tw)),
				labl=colnames(tw)
		)
		
		biomass.final<-bind.fv(x=biomass.final,
				y=data.frame("Biomass"=apply(res$v,1,mean,na.rm=TRUE)),
				desc=paste("Mean neighbourhood",funtype,"over types"),
				labl="MeanBiomass",
				preferred="Biomass"
		)		
		# a frequency weighted mean instead of just a mean, w=freqs/sum(freqs)
		#Iw=apply(res$v,2,weighted.mean,w=w,na.rm=TRUE), 
	}
	
	# if target type given add the values for the target type
	else
	{
		biomass.final<-bind.fv(x=biomass.final,
				y=data.frame("Biomass"=res$v[,1]),
				desc=paste(funtype,"around type",target),
				labl="Biomass",
				preferred="Biomass"
		)
	}
	
	# attach the frequencies too
	attr(biomass.final,"frequencies")<-freqs(X[res$included])
	
	# and some notes
	attr(biomass.final,"neighbourhoodType")<-res$ntype
	attr(biomass.final,"note")<-res$note
	
	# point values
    attr(biomass.final,"point.values")<-res$point.values2
	
	# return 
	biomass.final
}

###############################################################################

#eof
