
###############################################################################
# Spatial Simpson index
#
#
# Author: Tuomas Rajala <tarajala@maths.jyu.fi>
# last change: 060909
###############################################################################

simpsonF<-function(X, r=NULL, ...)
#Simpson index for graphs, with possibly various range-parameters
{
	# check that X is multitype ppp-object
	verifyclass(X, "ppp")
	if(length(levels(X$marks))<2) stop("Use only on a multitype point pattern (data.frame-marks not yet supported).")
	
	# the main calc function: returns the typewise mean of (deg_i(o)/deg)^2
	res<-segregationFun(X, r=r, fun="simpson", ...)
	
	# calc the non-spatial (global) value
	aspat<-simpson.index(X,spatial=FALSE)
	
	#TODO: the CSR values: possibly not right, the mean degree is not included
	theo<-rep(aspat, length(res$parvec))
	
	# calc the spatial value
	S <- 1 - unname(rowSums(res$v))
	
	# create the fv-object
	simpson.final<-fv(data.frame(theo=theo, par=res$parvec, S=S),
			argu="par",
			alim=range(res$parvec),
			ylab=substitute(S,NULL),
			desc=c("CSR values","Parameter values","Spatial Simpson index"),
			valu="S",
			fmla=".~par",
			unitname=res$unitname,
			fname="Spatial Simpson index"
	)
	
	
	
	# include also the typewise values of which the index is a summary
	attr(simpson.final,"typewise")<- -res$v
	
	# add the global index value too
	attr(simpson.final,"Aspatial Simpson index")<-aspat
	
	# add the note about neighbourhood definition
	attr(simpson.final,"note")<-res$note
	
	# return
	simpson.final
}

####################################################################################################
#Simpson index, just one value

simpson.index<-function(X, spatial=FALSE, ...)
{
	#the traditional aspatial Simpson index 1-D
	if(!spatial)
	{			
		sum0<-summary(X)
		ints<-sum0$marks[,3]
		m<-union(X$marks,NULL)
		pii<-ints/sum(ints)
		N<-X$n
		n<-pii*N		
		S<- 1 - sum( n*(n-1)  )/(N*(N-1))#Simpson index of diversity
		names(S)<-"Non-spatial Simpson index"
	}
	if(spatial)
	{			   #spatial Simpson index for a set of edgelists	
		S<-simpsonF(X, ...)$S
		names(S)<-"Spatial Simpson index"
	}
	S
}
