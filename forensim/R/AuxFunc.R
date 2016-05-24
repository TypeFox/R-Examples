#############################################
#Auxilary Functions#
#Hinda Haned,  december 10th 2008
#############################################

# returns the  position and the value of the maximum of a vector 

findmax<-function(vec)
{
	res0<-matrix(0,nrow=1,ncol=2)
	#res0<-data.frame(res0)
	if(all(vec==0))
	{
		xmax<-0
		maxval<-0
	}
	else{
	xmax<-which(vec==max(vec))
	maxval<-(vec[xmax])
	}
	#print(round(xmax))
	
	res0[1,1]<-xmax
	res0[1,2]<-signif(maxval,2)
	
	colnames(res0)<-c('max','maxval')
	rownames(res0)<-names(vec)
	return(res0)
}




#takes a matrix or a data frame in the format of the  Journal of Forenis Sciences for genetic data, and returns a list of markers with their allele frequencies,
#skips the NA markers
naomitab<-
function(tab)
{
	lem<-length(colnames(tab))-1 
	liste<-vector('list',lem)
	for(o in 2:(lem + 1))
	{
		tmp<-na.omit(tab[,c(1,o)])#data.frame
		rownames(tmp)<-NULL
		tmp1<-as.factor(tmp[,1])
		tmp2<-tmp[,2]
		names(tmp2)<-tmp1
		
		liste[[o-1]] <- tmp2
		
	}
	names(liste) <- colnames(tab)[-1]
	return(liste)
}



# finds the corresponding allele frequencies  of a simumix object from a tabfreq object

findfreq <-
function(mix,freq, refpop=NULL)
{
	pop <- freq@pop.names
	mix.all <- mix@mix.all
	mixmark <- mix@which.loc
	mark2 <- freq@which.loc
	
	
	
	#freq contains only one population, wuthout a name
	
	 if(!all(mixmark %in% mark2))
	 {	
	 stop('given objects for mix and freq are not compatible, freq must be the object from which the mixture was generated')
	 }
	#if freq contains  only one population, withe null name
	if(is.null(pop))
	{
		freq <- freq@tab
		res <- vector('list',length(mixmark))
		names(res) <- mixmark
		
		for( j in mixmark)
		{
			A <- mix.all[[j]]
			
			tmp0 <- freq[[j]][A]
			names(tmp0) <- A
			y0 <- which(is.na(tmp0))
			tmp0 <- replace(tmp0, y0, 0)
			res[[j]] <- tmp0
							
				
		}
			
		return(res)
	}
	
	#if freq contains multiples populations (or one pop with a name), either one is chosen, specified with refpop, or all populations in the mixture are considered
	if(!is.null(pop))
	{
		freqpop <- freq@tab
		#print(refpop)
		#no reference pop specified
		if(length(pop)==1)
		{
			respop <- vector('list', 1)
			names(respop) <- pop
			freq2 <- freqpop[[pop]][mixmark]
			res2 <- vector('list',length(mixmark))
			names(res2) <- mixmark
					
			for(k in mixmark)
			{	#special case when alleles in the mixture are not found un the reference population
				A2 <- mix.all[[k]]
				tmp <-freq2[[k]][A2]
					
				names(tmp) <- A2
				y <- which(is.na(tmp))
				tmp <- replace(tmp,y, 0)
				res2[[k]] <- tmp
			}
				
			respop[[pop]] <- res2
		}
		
	
		#no need to give a name for refpop
		
		if(length(pop)>1)
		{
			if(is.null(refpop))
			{	
				#populations in freq and mix mus then be the same
				stop("no reference population was given")
			}
				
			if(!is.null(refpop))
			{
				respop <- vector('list', 1)
				names(respop) <- refpop
				if(!(refpop %in% pop))
				{
					stop('uknown reference population')
				}
					
				freq2 <- freqpop[[refpop]][mixmark]
				res2 <- vector('list',length(mixmark))
				names(res2) <- mixmark
					
				for(k in mixmark)
				{	#special case when alleles in the mixture are not found un the reference population
					A2 <- mix.all[[k]]
					tmp <-freq2[[k]][A2]
					
					names(tmp) <- A2
					y <- which(is.na(tmp))
					tmp <- replace(tmp,y, 0)
					res2[[k]] <- tmp
					
					
				}
				
				respop[[refpop]] <- res2
			}
		}
		return(respop)
	}
}





# the number of possible combinations of   r objects  among c with repetitions
Cmn <- 
function(m,n)
{
	factorial(n+m-1)/(factorial(n-1)*factorial(m))
}

#gives all the possible permutations of the  Cmn(m,n)  combinations , calls a n function recurs 
comb <- function(m, n){
nL<-round(Cmn(m,n))
matR <- matrix(0,nrow=nL,ncol=n)
recurs<-function(m,n,matR,nbLigne,nbCol)
{
.C('recurs',as.integer(m),as.integer(n),matR=as.integer(matR),as.integer(nbLigne),as.integer(nbCol), PACKAGE="forensim")
}
loc<-recurs(m,n,matR,nL,n)
matrix(loc$matR,ncol=n)
}

	
	
	
	
	
# number of alleles in a mixture per locus, or among loci  (default)
nball <-
function(mix, byloc=FALSE)
{
	if(byloc)
	{
		
		return(sapply(mix@mix.all, length))
		
	}
	
	else
	{
		return(sum(sapply(mix@mix.all, length)))
	}

}





