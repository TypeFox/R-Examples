# Hinda Haned, December 2008
# haned@biomserv.univ-lyon1.fr

# simulate  allele frequencies from  a Dirichlet distribution
# results are given in the format of the Journal of Forensic Science for genetic data 

simufreqD <- 
function(nloc=1,nal=2, alpha=1)
{
	#from the  gtools library : function used to genrate allele frequencies in a population,
	rdirichlet <- function(n, a)
	{
        le <- length(a)
        x <- matrix(rgamma(le * n, a), ncol = le, byrow = TRUE)
        sm <- x %*% rep(1, le)
        x/as.vector(sm)
    }
	
	if(length(nal)!=1 & length(nal)!=nloc)
	{
		stop("nal must be of length ", nloc)
	}
	
	
	
		
	maxna <- max(nal)
	p <- matrix(0, nrow = maxna, ncol= nloc)
	
	#if alpha is a vector :
	if(is.vector(alpha))
	{
		# it can either be an integer specifying the alpha parmeter for all alleles for all markers
		if(length(alpha)==1)
		{	
			if(length(nal)==1)
			{
				nal   <- rep(nal,nloc) 
				
			}
			for(i in 1:nloc)
			{
						#tmp <-rdirichlet(1, rep(alpha[i], nal[i]))
						
						if(nal[i]==maxna)
						{
							p[,i] <- rdirichlet(1, rep(alpha, nal[i]))
						}	
						
						else
						{	#if the number of alleles is different per marker
							k <- maxna-nal[i]
							p[,i] <- c(rdirichlet(1, rep(alpha, nal[i])),rep(NA,k))
						}	
			}
		}
		
		#or it can be a vetor giving for alle markers, the alpha parameters
		if(length(alpha)!=1)
		{
			#nal can be of length one with a vector for alpha, but if both alpha and nal are  vectors, they must have the same length
			
			s<- sum(nal)
			# if(identical(length(nal),1) & !(identical(length(alpha),sum(nal))))
			if(length(nal)==1 & length(alpha)!=s)
			{
				
				stop("alpha must be of length ", sum(nal))
			}
			if(length(nal)==1)
			{
				nal   <- rep(nal,nloc) 
				
			}
			
		
			for(i in 1:nloc)
			{
					#tmp <-rdirichlet(1, rep(alpha[i], nal[i]))
				if(nal[i]==maxna)
				{
					#p[,i] <<- rdirichlet(1, rep(alpha, nal[i]))
					
					p[,i] <- rdirichlet(1, alpha)
				}	
				else
				{	#if the number of alleles is different per marker
					k <- maxna-nal[i]
					p[,i] <- c(rdirichlet(1, alpha),rep(NA,k))
				}	
			}
		}
	}
	
	if(is.matrix(alpha))
	{	#a matrix must specify the variance parmeters for each marker, wether all markers share the same number of alleles or not (the two cases are treated)
		#if the markers don't share thes same number of alleles , NA's must be  intoduced in the  parameters  matrix
		if(length(nal)==1)
		{
			nal   <- rep(nal,nloc) 
		}
		if(all(dim(alpha)==c(nloc,maxna)))
		{
		
			for(i in 1:nloc)
			{
					#tmp <-rdirichlet(1, rep(alpha[i], nal[i]))
				if(nal[i]==maxna)
				{
					#p[,i] <<- rdirichlet(1, rep(alpha, nal[i]))
					
					p[,i] <- rdirichlet(1, alpha[i,])
				}	
				else
				{	#if the number of alleles is different per marker
					k <- maxna-nal[i]
					v <- na.omit(alpha[i,])
					#print(v)
					p[,i] <- c(rdirichlet(1, v),rep(NA,k))
				}	
			}
		}
		
		
		else
		{
			stop("alpha must be a ",nloc,"x", maxna, " matrix")
		
		}
	}
		
	
	p <- cbind(as.factor(1:maxna), p)		
	colnames(p) <- c('Allele', paste('Marker',1:nloc,sep=''))
	rownames(p) <- 1:maxna
	return(signif(p,2))
}