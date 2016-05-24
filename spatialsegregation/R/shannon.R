# Spatial Shannon index functional
#
#
# Author: Tuomas Rajala <tarajala@maths.jyu.fi>
###############################################################################

shannonF<-function(X, r=NULL, v2=FALSE, ...)
#Shannon index for graphs, with possibly various range-parameters
{
	# check that X is multitype ppp-object
	verifyclass(X, "ppp")
	if(length(levels(X$marks))<2) stop("Use only on a multitype point pattern (data.frame-marks not yet supported).")
	
	# the main calc function: it calculates only the pi_tau-vector
	res<-segregationFun(X=X, r=r, fun="shannon", funpars=ifelse(v2,1,0), ...) 
	
	# calc the aspatial (global) entropy i.e. shannon index
	eglobal <- -shannon.index(X,spatial=FALSE)
	
	# a function to calculate  -1*sum( pi_tau*log(pi_tau) )
	f<-function(pvec) 
	{
		E1<-0
		pvec<-pvec/sum(pvec)
		S<-length(pvec)
		ok<-pvec>0
		E1 <- sum(pvec[ok]*log(pvec[ok],base=S))
		(1-E1/eglobal)
	}
	
	# if we take the log base as the individual degree instead of total S -
	#   then the result is the mean
	if(v2)
	{
		H<-apply(-res$v,1,mean)
		desc<-"Spatial Shannon index, v2 (check documentation)"
	}
	# the base is S: calc the E(o)/E
	else
	{
		H<-unname(apply(res$v,1,f))
		desc<-"Spatial Shannon index"
	}
		
	# the CSR values
	theo<-rep(0, length(res$parvec))
	
	# create the fv-object
	shannon.final<-fv(data.frame(theo=theo, par=res$parvec, H=H),
			          argu="par",
					  alim=range(res$parvec),
					  ylab=substitute(H,NULL),
					  desc=c("CSR values","Parameter values",desc),
					  valu="H",
					  fmla=".~par",
			      unitname=res$unitname,
					 fname=desc
					 )
	
	# add the typewise pii_tau values of which the index is a summary
	attr(shannon.final,"typewise")<- res$v
	
	# add the global index value too
	attr(shannon.final,"Aspatial Shannon index")<--eglobal
	
	# add the note about neighbourhood definition
	attr(shannon.final,"note")<-res$note
	
	# return
	shannon.final
}


####################################################################################################

#Shannon index

shannon.index<-function(X, spatial=FALSE, ...)
{
	#the traditional aspatial Shannon index
	if(!spatial)
	{  
		sum0<-summary(X)
		ints<-sum0$marks[,3]
		m<-union(X$marks,NULL)
		pii<-ints/sum(ints)
		H<- (-1)*sum( pii * log(pii, length(m)),na.rm=TRUE)
		names(H)<-"Non-spatial Shannon index"
	}
	if(spatial)
	{			   #spatial Simpson index for a set of edgelists	
		H<-shannonF(X=X, ...)$H
		names(H)<-"Spatial Shannon index"
	}
	H
}
