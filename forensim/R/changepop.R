# Hinda Haned, April 2009
# haned@biomserv.univ-lyon1.fr

changepop <- 
function(obj,oldpop,newpop)
{
	if(is.simugeno(obj) || is.tabfreq(obj))
	{
		pop <-obj$pop.names
	
	
		if(length(newpop)!=length(oldpop))
		{
			stop('newopop and ooldpop must have the same length')
		}
		if(!all(oldpop %in% pop)) stop('uknown population names in oldpop')
		
		
		if(length(pop)==1 & length(newpop)==1) #and le(oldpop)==1
		{	#only one pop
			if(is.tabfreq(obj))
			{
				#print('ici')
				obj$pop.names <- factor(newpop)
				names(obj$tab) <- newpop
			}
			if(is.simugeno(obj))
			{
				obj$pop.names <- factor(newpop)
				names(obj$tab.freq) <-newpop
				obj$popind <- as.factor(rep(newpop,obj$nind)) #as all ind belong to the same population
			}
					
		
					
		}
		
		
		if( length(pop) > 1 & length(newpop) >=1)# and length(oldpop) >1
		{	#more than one pop nae to change
			if(is.tabfreq(obj))
			{
				tmp <- as.character(obj$pop.names)
				index <- which(oldpop %in% tmp)
				tmp[index] <- newpop #2 or more changes
				obj$pop.names <- as.factor(tmp)
				names(obj$tab) <- tmp
			}
			
			if(is.simugeno(obj))
			{
				tmp <- as.character(obj$pop.names)
				index <- which(oldpop %in% tmp)
				tmp[index] <- newpop #2 or more changes
				obj$pop.names <- as.factor(tmp)
				names(obj$tab.freq) <- tmp
			
			}
			
			
		}
	}

	if(is.simumix(obj))
	{
		tmp1 <- table(obj$popinfo)[oldpop]
		#print(tmp1)
		names(tmp1) <- NULL
		#print('tmp1'); print(tmp1)
		tmp2 <- as.character(obj$popinfo)
		index2 <- which(oldpop==tmp2)
		tmp2[index2] <- rep(newpop,tmp1)
		obj$popinfo <- as.factor(tmp2)
	}	
	
	return(obj)	#if(is.simumix(obj))
	
	
}