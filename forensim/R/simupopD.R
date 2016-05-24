# Hinda Haned, December 2008
# haned@biomserv.univ-lyon1.fr

#simulate  populations allele frequencies from a Dirichlet Distribution


simupopD <-
function( npop=1, nloc=1, na = 2,  globalfreq = NULL, which.loc=NULL,  alpha1, alpha2 = 1) 
{
    
	# globalfreq is  in JFS format
	#from the  gtools library : rdirichlet
	rdirichlet <- function(n, a)
	{
        le <- length(a)
        x <- matrix(rgamma(le * n, a), ncol = le, byrow = TRUE)
        sm <- x %*% rep(1, le)
        x/as.vector(sm)
    }
	
	
	
	if(identical(npop,0) || is.null(npop))
	{
		stop("invalid number of populations")
	}	
	if(identical(nloc,0) || is.null(nloc))
	{
		stop("invalid number of markers")
	}
	if(length(na)!=1 & length(na)!=nloc)
	{
			stop("na must be of length ", nloc)
	}
	
	if(length(alpha1)!=npop)
	{
		stop("npop and alpha1 must have the same length ")
	}
	
	#default is the uniform distribution for allele frequencies  in the globalfreqal population
if(is.null(globalfreq))
	{
		glob <- as.matrix(simufreqD(nloc,na,alpha2))
		loki<-colnames(glob)[-1]
		globF<-tabfreq(glob) #for the final result 
		tab <- vector('list',npop)
		mat <- matrix(0,nrow = max(na),ncol = nloc)
		colnames(mat) <- paste('Marker', 1:nloc,sep="")
		rownames(mat) <- 1:max(na)
		#for each pop
		for(j in 1:npop)
		{
			#each variance parameter correspond to a subpopulation
				for(i in loki)
				{	#for each marker 
					tmp <- na.omit(glob[, i])  
					varp <- tmp * (1 - alpha1[j])/alpha1[j]  #variance parameter
					if(length(tmp)==max(na))
					{				
						temp1 <- rdirichlet(1, varp)
					}
					else
					{
						k <- max(na)-length(tmp)
						temp1 <- c(rdirichlet(1,varp),rep(NA,k))
						#temp1: AF for the ieme locus in the jeme pop

					}
					
					mat[,i]<-signif(temp1,4)
				}	
			Allele <- rownames(mat)
			
			mat2 <- cbind.data.frame(Allele,mat)
			
			tab[[j]] <- mat2
		}
		pop.names <- as.factor(paste("pop",1:npop, sep=""))
		res <- tabfreq(tab,pop.names)
		   
   }
   
   else
   {
		tab <- vector('list',npop)
		Allele <- globalfreq$Allele	
		loc <- vector('list',npop)
		if(is.null(which.loc))
		{
			#all markers in the ref database are considered 
			globF <- tabfreq(globalfreq)
			k <- globF@which.loc
		}
		
		else
		{
			k <- which.loc
			foo <- globalfreq[,c('Allele',k)]
			globF <- tabfreq(foo)
		}
		
		
		
		
		for(i in 1:npop)
		{
			listemp <- vector('list', length(k))
			names(listemp) <- k
			
			for(j in k)
			{					
				tmp <- globF@tab[[j]] 
				locA <- names(tmp)
				#print(locA)
				varp <- tmp * (1 - alpha1[i])/alpha1[i]  #variance parameter
				temp1 <- rdirichlet(1, varp)
				temp1 <- as.vector(temp1)
				names(temp1)<-locA
				listemp[[j]]<-signif(unlist(temp1),4)
			}
			#print(listemp)
			tab[[i]]<- listemp
		}
		names(tab)<- paste("pop",1:npop, sep="")
		res <- new('tabfreq')
		res@tab <- tab
		res@pop.names <- as.factor(paste("pop",1:npop, sep=""))
		res@which.loc<- globF@which.loc
		
	}
   
	final <-vector('list',2)
	final[[1]]<-globF
	final[[2]]<- res
	names(final) <-c("globfreq", "popfreq")
	return(final)
   
}

