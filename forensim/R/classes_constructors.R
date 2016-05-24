# Hinda Haned, December 2008
# haned@biomserv.univ-lyon1.fr
##############################################
# forensim  classes constructors definitions 
##############################################
#making sure that the simulations are reproducible
#simugeno class constructor :  simugeno function 
simugeno <-
function(tab,which.loc=NULL,n=1)
{
	if(!is.tabfreq(tab))
	{
		stop("tab must be a tabfreq object")
	}
	X <- tab@tab
	M <- tab@which.loc
	#cheking that marker.names are correct
	if(is.null(which.loc))
	{	
		which.loc <- M
	}
	else
	{	
		if(identical(which.loc, character(0)))
		{
			stop("invalid markers names")
		}
		
		if(!all(which.loc %in% M))
		{	
			stop("uknown markers  chosen, please check the markers names")
		}
	}
	popfac <- tab@pop.names
	#n <- as.integer(n)
	#locn is a local variable defining the number of genotypes 
	if(is.null(n))
	{
		locn <- 1
	}
	
	if(any(n  < 0) )
	{
		stop('n must be a positive integer, or a list of positive integers')
	}	
	if(identical(length(n),1) & identical(n,0))
	{
		stop("0 is not a valid value for n ")
	}
	
	# one population 
	if(is.null(popfac) & length(n) >= 1 )
	{
		popind <- NULL
		if(length(n) > 1) 
		{
			locn <- n[1]
			if(locn == 0)
			{
				stop("0 is not a valid value for n ")
		
			}
				warning("only the first element of n will be considered")
				
		}	
		if(length(n)==1  )
		{
			locn <- n 
		}
		
		if(is.null(n))
		{
			locn <- 1
		}
		
		geno <- matrix(0,ncol=length(which.loc),nrow=locn)
		rownames(geno) <- paste("Ind",1:locn,sep='')
		colnames(geno) <- which.loc
				
		for(i in which.loc)
		{
			a1 <- sample(names(X[[i]]),locn,prob=X[[i]],replace=TRUE)
			a2 <- sample(names(X[[i]]),locn,prob=X[[i]],replace=TRUE)
			geno[,i] <- paste(a1,a2,sep="/")
		}
		
		
		freq <- X[which.loc]#if one pop
	}
	
	# if multiple pops
	if(!is.null(popfac))
	{
		#verifications
		if(length(n) > length(popfac))
		{
			locn <- n[1:length(popfac)]
			if(all(n==0))
			{
				stop("n is empty for all populations ")
			}
			warning("n contains too many elements, only the first ",length(popfac), " are considered" )
		}
		
		if(length(n) < length(popfac))
		{	
			stop("n must be of length ",length(popfac))
		}
		#verif end
		if(identical(length(n),length(popfac)))
		{
			locn <- n
		}
		names(locn)<-popfac
		#print('popfac');print(popfac)
		#print(length(popfac))
		freq <- vector('list',length(popfac))
		names(freq)<- popfac
		nombL <- length(which.loc)
		names(locn) <- popfac
		X <- tab@tab
		geno0 <- NULL
		geno <- NULL
		popind <- NULL
		tempfreq <- NULL
		if(length(popfac)==1)
		{
			popind<-rep(popfac,n)
			for(i in popfac)
			{
				geno0 <- matrix(0,ncol=nombL,nrow=locn[i])#changes at each iteration, different number of sampeled individuals
				for(m in which.loc)
				{	
					colnames(geno0) <- which.loc
					if(locn[i]==0)
					{
						geno0[ , m] <- 0
					}
					else
					{
						avec1 <- sample(names(X[[i]][[m]]),locn[i],prob=X[[i]][[m]],replace=TRUE)
						avec2 <- sample(names(X[[i]][[m]]),locn[i],prob=X[[i]][[m]],replace=TRUE)
						ty <- paste(avec1,avec2,sep='/')
						geno0[ , m] <- ty 
						tempfreq <- c(tempfreq,X[[i]][m])#concatenation for all markers for one pop
					}
				}
					
				geno <- rbind(geno,geno0)
				freq[[i]]<- tempfreq
				tempfreq <- NULL#free 
				#names(freq[i]) <- i
				#print(freq[i])
			}
		}
		#popind
		if(length(popfac) > 1)
		{
			for(i in popfac)
			{
				geno0 <- matrix(0,ncol=nombL,nrow=locn[i])#changes at each iteration, different number of sampeled individuals
				#print(rep(levels(popfac)[i],locn[i]))
				popind <- c(popind,rep(i,locn[i]))
				#print(popind)
				for(m in which.loc)
				{	
					colnames(geno0) <- which.loc
					if(locn[i]==0)
					{
						geno0[ , m] <- 0
					}
						
					else
					{
						avec1 <- sample(names(X[[i]][[m]]),locn[i],prob=X[[i]][[m]],replace=TRUE)
						avec2 <- sample(names(X[[i]][[m]]),locn[i],prob=X[[i]][[m]],replace=TRUE)
						ty <- paste(avec1,avec2,sep='/')
						geno0[ , m] <- ty 
						tempfreq <- c(tempfreq,X[[i]][m])#concatenation for all markers for one pop
					}
				}
						
					#freq[[g]] <<- tempfreq
					geno <- rbind(geno,geno0)
					freq [[i]]<- tempfreq
					tempfreq <- NULL#free 
			}
		}

	popind <- as.factor(popind)
	}
	ID <- paste('ind',1:sum(locn),sep='')
	rownames(geno) <- ID
	res <- new('simugeno')
	res@tab.freq <- freq 
	res@nind <- sum(locn) 
	res@which.loc <- which.loc
	res@tab.geno <- geno
	res@pop.names <- popfac
	res@popind <- popind 
	res@indID <- ID
	return(res)
}





# simumix class constructor: simumix function 
simumix <- 
function(tab,which.loc=NULL,ncontri=1)
{

	if(!is.simugeno(tab))
	{
		stop("\n tab must be a simugeno object\n")
	}
	
	
	M <- tab@which.loc
	
	# cheking marker names
	if(is.null(which.loc)) 
	{
		which.loc <- M
	}
	
	#if not null names, than check for possible errors
	else
	{
		if(identical(which.loc, character(0)))
		{
			stop("markers names must be a character string")
		}
		
		if(!all(which.loc %in% M))
		{	
			#if the user's markers are not the same that those of the tabfreq object 
			stop("uknown marker names")
		}
	}
	
	
	popfac <- tab@pop.names
	popind <- tab@popind
	nsimugeno <- tab@nind
	prof <- tab@tab.geno
	tabfreq <- tab@tab.freq
	nmark <- length(which.loc)
	mix.all <- vector('list',nmark)
	names(mix.all) <- which.loc
	indID <- tab@indID
	N <- ncontri 
		
	if(sum(ncontri) > nsimugeno)
	{
		stop("too much contributors in the mixture")
	}	
	if(any(ncontri < 0))
	{
		stop("ncontri must be a positive integer, or a list of positive integers")
	}
	if(is.null(popfac))
	{
			if(identical(ncontri, 0))
			{
				stop("0 is not a valid value for ncontri ")
			}
			if(length(ncontri) > 1)
			{
				
				N <- ncontri[1]
				if(N==0)
				{
					stop("0 is not a valid value for ncontri ")
		
				}
				warning("only the first element of ncontri will be considered")
				
			}	
			
			
			#sampling individuals
			
			samp0 <- sample(indID,N,replace=FALSE)#chosen profiles to ssimulate the mixture
			mixprof <- prof[samp0,which.loc] 
			mixprof <- as.matrix(mixprof)
			indnames  <- samp0
			if(length(which.loc)>1 & N==1)
			{
				mixprof <- t(mixprof)
				#colnames(mixprof) <- which.loc	
			}
			rownames(mixprof) <- indnames
			popinfo <- NULL
			colnames(mixprof) <- which.loc	
			#print(mixprof)
			
	}
		
	#multiple populations case
	if(!is.null(popfac))
	{	
		
		if(length(ncontri) <= 1)
		{
			if(!identical(length(popfac),length(ncontri)))
				stop("ncontri must be of length ", length(popfac))
		}
			
		if(length(ncontri) > 1 & all(ncontri==0))
			stop("ncontri is empty for all populations ")
		
		if(length(ncontri) > length(popfac))
		{
			warning("only the ", length(popfac), " elements of ncontri are considered")
			ncontri <- ncontri[1:length(popfac)]
		}
	
		if(length(ncontri) < length(popfac)){
			stop("ncontri must be of length ",length(popfac))}
		
		#case where the number of conributors is greater of the profiles present in the simugeno object
		names(N) <- popfac
		#print(N)
		temprof <- NULL
		temp0 <- NULL
		popinfo <- NULL
		poptemp <- as.character(unique(popind))
		#print(poptemp)
		for(p in popfac)
		{
			#print(N[p])
			if(N[p]!=0)
			{
				samp1 <- sample(indID[popind == p  ],N[p],replace=FALSE)#chosen profiles to ssimulate the mixture
				temp0 <- c(temp0,samp1)
				temp1 <- prof[samp1,]
				temprof <- rbind(temprof,temp1) 
				popinfo <- c(popinfo,rep(names(N[p]), N[p]))
			}				
					
		}
	indnames <- temp0
	mixprof  <- temprof
	mixprof <- mixprof[,which.loc]
	mixprof <- as.matrix(mixprof)
	if(ncol(mixprof)==1)
	{
		mixprof<-t(mixprof)
		rownames(mixprof)<- indnames
		colnames(mixprof)<-which.loc
	}
	else{
	rownames(mixprof) <- indnames
	colnames(mixprof) <- which.loc
	}
	popinfo <- as.factor(popinfo)
	#print('apres')
	}
		
	# alleles present in the mixture
	 for(i in which.loc)
	 {
			 loc <- sort(unique(unlist(strsplit(mixprof[,i],'/'))))
			 temp.loc<- as.character(sort(as.numeric(loc)))
			 mix.all[[i]] <- temp.loc
	 }
	
	res <- new('simumix')
	
	res@mix.all <- mix.all
	res@mix.prof <- mixprof
	res@ncontri <- sum(N)
	res@which.loc <- which.loc
	res@popinfo <- popinfo
	return(res)
}
	




#tabfreq class constructor: tabfreq function 
tabfreq <-
function(tab,pop.names=NULL)
{	
	
	#  single  population  case
	if(class(tab)!= "list" )#if it's not a list, it must be a data.frame or a matrix
	{	
		res <- vector('list',1)
		if(class(tab) != "data.frame" & class(tab)!="matrix")
		{#if tab is not a matrix or data.frame 
			stop("tab must be a of types  matrix or  data.frame")
		}
		
		if (is.null(colnames(tab))) 
		{
			stop("tab columns have no name")
		}	
		
		#cheking pop.names : one cold want to specify the population name, even for a single population
		if(length(pop.names)  > 1)
		{	
			if(!is.factor(pop.names))
			{
				stop("a factor is expected for pop.names")
			}
			warning("only the first element of pop.names is considered")
			popfac <- pop.names[1]
			names(res) <- popfac
			res[[popfac]] <- naomitab(tab)
			
		}
		
		
		if(length(pop.names) == 1)
		{	
			
			if(!is.factor(pop.names))
			{
				stop("a factor is expected for pop.names")
			}
			
			popfac <- pop.names
			names(res) <- popfac
			res[[popfac]] <- naomitab(tab)
			
		}
		
		
		if(is.null(pop.names))
		{
			popfac <- pop.names
			res <-0 
			res <- (naomitab(tab))
			
		}
		mark<- colnames(tab)[-1]
		
	}
	
	
	
	# multiple populations  case 
	if(class(tab) == "list")
	{	
		p <- length(tab)
		res <- vector('list',p)
		mm <- vector('list',p)
		popfac <- pop.names
		if(!is.factor(pop.names)){
			
			stop("a factor is expected for pop.names")
		}	

		
		if(length(pop.names) != p || length(pop.names) < p)
		{
			stop("tab and pop.names must have the same length")
		}
		
		if(length(pop.names) > p){
			popfac <- pop.names[1:p]
			warning("only the first ", p," elements of pop.names are considered")
		}
		
		
		for(i in popfac)
		{
			if(class(tab[[i]]) != "data.frame" & class(tab[[i]]) != "matrix")
			{
				stop("all elements of tab must be of types matrix or data.frame")
			}
			
			if (is.null(colnames(tab[[i]]))) 
			{
				stop("columns of element ", i, " of tab have  no name ")
			}
					
			res[[i]] <- naomitab(tab[[i]])
			mm[[i]] <- colnames(tab[[i]])[-1]
		}
		
		names(res) <- popfac
		
		#checking markers between populations
		maxloc <- max(sapply(mm,function(y) length(y)))
		lpop <- length(popfac)
		comloc <- names(which(table(unlist(mm))==lpop))
		if(length(comloc)==0)
		{
			stop("There is no shared markers between the ", lpop," populations. Please check your data.")
		}
	
		if(length(comloc) != 0 )
		{
			if(length(comloc) < maxloc)
			{
				warning("Only ", length(comloc)," markers are common between the ", lpop," populations. You might want to check your data.")
				for(w in popfac)
				{
					res[[w]]<-res[[w]][comloc]
					names(res[[w]]) <- popfac
				}
				
				mark <- comloc
			}
			if(length(comloc)==maxloc)
			{
				mark <- mm[[1]]#for example, no need to complicate
				#if condition not met "res"  doesn't change
			}
		}
		
	}
	
	#object tabfreq to construct
	freq <- new('tabfreq')
	freq@pop.names <- popfac
	freq@tab <- res
	freq@which.loc <- mark
	#freq@theta <- theta
	return(freq)
}



# Alias Method
as.simugeno <- simugeno
as.simumix <- simumix
as.tabfreq <- tabfreq
