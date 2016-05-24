# Hinda Haned, April 2009
# haned@biomserv.univ-lyon1.fr

PE <-
function(mix, freq, refpop=NULL, theta=0,byloc=FALSE)
{
	
	popinfo <-unique(mix@popinfo)
	freqinfo <- freq$pop.names
	loc <- mix@which.loc
	pe.loc<-matrix(0,nrow=1,ncol=length(loc))
	colnames(pe.loc) <- loc
	rownames(pe.loc) <- "PE_l"
	if(is.null(popinfo))
	{
		af <- findfreq(mix,freq)
		#print(af)
	}
	
	if(!is.null(popinfo))
	{	
		if(length(freqinfo)==1)
		{
			af <- findfreq(mix,freq)[[popinfo]]
		}
		
		if(length(freqinfo) >1)
		{
			#
			af <- findfreq(mix,freq, refpop)[[refpop]]
		}
	}	#if(length(popinfo>1))
	
	
	if(theta==0)
	{
		for( l in  loc)
		{
			#print(af[[l]])
			pe.loc[1,l]<-1-sum(af[[l]])^2
		}
		
		
		if(byloc)
		{
			return(signif(t(pe.loc),4))
		}
		
		else
		{
			PE <- 1-prod(1-pe.loc)
			names(PE) <- "PE"
			return(signif(PE,6))
		}
	}
	if(theta!=0)
	{
		if (theta >= 1 || theta < 0) 
		{
			stop("'theta' must be a number in [0,1[" )
		}
		
		for( l in  loc)
		{
			pe.loc[1,l]<-1-sum(af[[l]])^2-theta*sum(af[[l]])*(1-sum(af[[l]]))
		}
		
		
		if(byloc)
		{
			return(t(pe.loc))
		}
		
		else
		{
			PE <- 1-prod(1-pe.loc)
			names(PE) <- "PE"
			return(signif(PE,6))
		}
	}
}








