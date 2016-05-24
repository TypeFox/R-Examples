# ISAR or local species richness summary for multitype spatial point pattern
#
# Desc: 
#   Estimate the number of types present in the neighbourhoods of points
#
#
# Neighbourhood definitions (parameter): 
#   Geometric (r)
#   k-nearest neighbours (k)
#	Gabriel (none)
#   Delauney triangulation (none)
#
#
#
# Author: Tuomas A. Rajala <tuomas.rajala@iki.fi>
#
#
# Last update: 201010
###############################################################################

isarF<-function(X, r=NULL, target=NULL, v2=FALSE, v3=FALSE, v4=FALSE, ... )
{
	# check that X is ppp-object
	verifyclass(X, "ppp")
	if(length(levels(X$marks))<2) stop("Use only on a multitype point pattern (data.frame-marks not yet supported).")
	# if no target given, calculate for all types
	if(is.null(target))
	{
		targeti <- 0
		valu   <- "ISARmean"
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
	
	#v2 logical if a degree weighted version should be calculated
	if(v2) funtype <- "Neighbour-count-weighted-ISAR"
	if(v3){v2<-2; funtype <- "Mass-weighted ISAR"; }
	if(v4){v2<-3; funtype <- "eISAR"; }
	else funtype <- "ISAR"
	
		
	# use the main calc function
	res<-segregationFun(X=X, fun="isar", r, funpars=c(targeti, as.integer(v2)), ...)
	
	# theoretical values in CSR: depends on the neighbourhood type
	ntype<-res$ntype
	mdeg<-function(l,k)c( pi*l*k^2, k, 4, 6)[charmatch(ntype, kGraphs)]
	
	# get the intensities
	sum0<-summary(X)
	omarks<-order(union(X$marks,NULL)) # right order of marks
	l<-sum0$marks[,3][omarks]
	#    calc the theoretical values, also for the degree weighted version
	theo<-NULL
	for(para in res$parvec)theo<-
				c(theo, 
                 sum(1-exp(-mdeg(sum(l),para)*l/sum(l) )) / ifelse(v2,mdeg(sum(l), para),1))
										
		
	# create the fv-object
	isar.final<-fv(data.frame(theo=theo,par=res$parvec), 
			       argu="par",
				   alim=range(res$parvec),
				   ylab=substitute(ISAR, NULL),
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
		colnames(tw)<-union(marks(X[res$included]),NULL)
		isar.final<-bind.fv(x=isar.final,
						    y=tw,
					  	 desc=paste("Typewise",funtype,"for type",colnames(tw)),
					     labl=colnames(tw)
					         )
					
		isar.final<-bind.fv(x=isar.final,
				            y=data.frame("ISARmean"=apply(res$v,1,mean,na.rm=TRUE)),
				         desc=paste("Mean",funtype,"over types"),
						 labl="ISARmean",
						 preferred="ISARmean"
				 	     	 )		
		# a frequency weighted mean instead of just a mean, w=freqs/sum(freqs)
		#Iw=apply(res$v,2,weighted.mean,w=w,na.rm=TRUE), 
	}
	
	# if target type given add the values for the target type
	else
	{
		isar.final<-bind.fv(x=isar.final,
				            y=data.frame("ISAR"=res$v[,1]),
						 desc=paste(funtype,"for type",target),
						 labl="ISAR",
					preferred="ISAR"
								)
	}

	# attach the frequencies too
	attr(isar.final,"frequencies")<-freqs(X[res$included])
	
	# and some notes
	attr(isar.final,"neighbourhoodType")<-res$ntype
	attr(isar.final,"note")<-res$note
	# add pointwise values
	attr(isar.final,"point.values")<-res$point.values
	# return 
	isar.final
}

###############################################################################
isar.index<-function(X, r=4, ntype="knn", ...)
{
	if(length(r)>1)stop("Use isarF for vector of parameter values.")
	I0<-isarF(X, r=r, ntype=ntype, ...)
	data.frame(meanISAR=I0$I, par=I0$par)
}