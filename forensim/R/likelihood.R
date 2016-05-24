#  January 9th  2009#
# Hinda Haned
#likelihood of the data, conditioning on x : the number of contributors to the mixture
# we follow here the notations of the general formula of likelihood rations interpretation in case of population substructure Curran et al , Journal of Forensic Sciences 1999
# x:  the number of contributors
#c: the number of distinct alleles in the mixture
#ui: unknown number of copies of allele Ai in the mixed stain
#r: the number of unconstrained alleles
#ri: the number of unconstrained alleles of type Ai 
#theta: Wright's Fst coefficient, correction for uncertainty on subpopulation allele frequencies in case of subdivision

#This function have several similarities with the function "Pevid.gen" of the forensic package, as it implemnets a particular 
#case of the latter

dataL <-
function(x=1,p, theta = 0 ) 
{
	
    if (!is.numeric(x) || is.na(x) || x < 0)
	{
        stop("'x' must be a positive integer")
    }
 
    if (!is.vector(p)) 
	{
        stop("a vector is expected for 'p'")
	}
    c <- length(p)
	#not enough contributors to explain the c alleles
    if(c > 2*x)
	{
		return(0)
	}
	
	else
	{
	
		if (round(sum(p)) > 1) 
		{
			stop("sum of allele frequencies must not exceed 1")
		}
		
		if (!is.numeric(p) || any(is.na(p)) || any(p < 0) || any(p > 1))
		{
			stop("alleles frequencies in  'p' must be numbers strictly between 0 and 1")
		}
	
	
   		if (theta >= 1 || theta < 0) 
		{
			stop("'theta' must be a number in [0,1[" )
		}
			
		r <- 2 * x - c
		
		
    
	    if (r>0)
		{ 
	        # all possible combinations of the unconstrained alleles  (r1, r2, ... , rc-1), where  
	        # sum(ri)=r and ri >= 0, ri <= r
	        ri <- comb(r, c)
		}
	    
		else
		{
	        ri <- matrix(rep(0, c), nrow = 1, byrow = TRUE)
		}
	    
	    # There are sum(ui) = 2*x  alleles in the mixture, not necessarily distinct
		#ui is the matrix giving all possible  numbers of copies of allele Ai in the stain
		#ui[,i]: number of alleles of type Ai, there are c possible types
	    ui <- matrix(0, nrow = nrow(ri), ncol = c)
	    ui <- ri + 1#+1 because one is already observed
		
	    if (x == 0)
		{
	        stop('The number of contributors x must be different from zero')
		}
		else
		{
		#if(2*x<c)
		    res1 <- factorial(2 * x)/prod((1 - theta) +  (0:(2 * x - 1)) * theta)
		} 
		
	    res2 <- rep(0, nrow(ui))
		
	    for(y in 1:nrow(ui))
		{ 
			#numerator
			num <- rep(0, c)
			for (i in 1:c)
			{
				if (ui[y, i] == 0)
				{
				#if ui=0
					num[i] <- 1  #0!
				}
				else
				{
					num[i] <- prod((1 - theta) * p[i] + (0:(ui[y, i] -1)) * theta)
				}  
			}	
			res2[y] <- prod(num)/prod(factorial(ui[y, ]))#prod d'un produit pour le double produit II
	    }

    return(res1 * sum(res2))
	}
}

#returns  a matrix, giving the  likelihood of the alleles conditional on the number of contributors for each marker 
#lik.loc
lik.loc <-
function(x=1, mix,freq,refpop=NULL,theta=NULL, loc=NULL)
{ 
	# call to Auxiliary Function findfreq()
	
	if(!is.null(loc))
	{
		mark <- loc
	}
	
	else
	{
		mark <- mix@which.loc
	}
	
	mark2 <- freq@which.loc
	mix.all <- mix@mix.all
	popfac <- freq@pop.names
	popmix <- unique(mix@popinfo)
	
	
	freqmix <- findfreq(mix,freq,refpop)
	# returns a list
	#single population case
	if(!all(mark %in% mark2))
	{
		stop('given objects for mix and freq are not compatible, freq must be the object from which the mixture was generated')
	}
	
	
	if(is.null(theta))
	{
		theta <-0
	}
	#if freq contains no populational  information
	if(is.null(popfac))
	{
			
			res <- sapply(mark, function(m) { 
			#print(m)
			#print(freqmix[[m]])
			dataL(x,freqmix[[m]], theta) })
		
	}
	#if freq contains only one population, no ambiguity		
	if(length(popfac)==1)
	{
		res<- sapply(mark, function(g) dataL(x,freqmix[[popfac]][[g]], theta))
	}
	
	if(length(popfac) > 1)
	{
		if(is.null(refpop))
		{
			stop('freq contains more than one population, please sepcify the reference population')
		}
			
		if(!is.null(refpop))
		{
		
			if(!(refpop %in% popfac))
			{
				stop('unknown population name')
			}
		
			#tmpfreq <- freq@tab						
			res <- sapply(mark, function(g) { dataL(x,freqmix[[refpop]][[g]], theta) })
			
		}
	}
	
	return(res)	
}

#Overall Likelihood, independence is assumed between markers (product of likelihoods)
lik <-
function(x=1, mix,freq,refpop=NULL,theta=NULL, loc=NULL)
{
	res <- lik.loc(x, mix,freq,refpop,theta, loc)
	return(prod(res))
}
