# January 9th  2009#
# Hinda Haned
#


#maximum likelihood estimation of the number of contributors, for each marker

likestim.loc<-
function(mix,freq,refpop=NULL,theta=NULL, loc=NULL)
{
	
	if(is.null(loc))
	{
			mark <- mix@which.loc
	}
	else
	{
			mark <- loc
	}
	#the maximum is searched in the discrete interval : 1: 6, more contributors is unlikely ?
	locres <- as.matrix(apply(sapply(1:6,function(i) lik.loc(i,mix,freq,refpop,theta, loc)),1,findmax))
	rownames(locres) <- c('max','maxval') 
	return(t(locres))
	
}
#likestm.loc(mix1,freq1)


#
# maximum likelihood , estimation of the number of contributors for all markers, independence is assumed between markers
#if not, uses the Balding et Nichols (1994)  correction for  population subdivision effect 

likestim<-
function( mix,freq,refpop=NULL,theta=NULL, loc=NULL)
{	
	tmp1 <- NULL
	for(h in 1:6)
	{
		tmp1 <- c(tmp1, lik(x=h, mix,freq,refpop,theta, loc))
	}
	return(findmax(tmp1))
	#return(tmp2)
}
#likestim((mix1,freq1)


