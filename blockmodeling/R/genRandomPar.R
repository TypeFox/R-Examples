"genRandomPar" <-
function(
k,#number of clusters/groups
n,#the number of units in each mode
seed=NULL,#the seed for random generation of partitions
mingr=1,	#minimal alowed group size
maxgr=Inf,	#maximal alowed group size
addParam = list(
	genPajekPar = TRUE, 	#Should the partitions be generated as in Pajek (the other options is completly random)
	probGenMech = NULL) #Here the probabilities for the 4 different mechanizems for specifying the partitions are set. It should be a numeric vector of length 4. If not set this is determined based on the previous parameter.
){
    if(is.null(addParam$probGenMech)){
    	if(is.null(addParam$genPajekPar)||addParam$genPajekPar) probGenMech <- c(1/3,1/3,1/3,0) else probGenMech <- c(0,0,0,1)
    } else probGenMech<-addParam$probGenMech
    if(!is.null(seed))set.seed(seed)
    nmode <- length(k)
    ver<-sample(1:4,size=1,prob=probGenMech)
      if(nmode==1){
        find.new.par<-TRUE
	while(find.new.par){
	  if(ver!=4){
		temppar<-integer(n)

		if(ver==1){
			temppar<-1:n%%k+1	
		}

		if(ver==2){
			temppar[1:k]<-1:k
			temppar[(k+1):n]<-k
		}

		if(ver==3){
			temppar[1:k]<-1:k
			temppar[(k+1):n]<-1+trunc(k*runif(n-k))
		}

		for(ii in n:2){
			jj<-trunc(ii*runif(1))
			temppar[c(ii,jj)]<-temppar[c(jj,ii)]
		}
          }else temppar<-sample(1:k,n,replace=TRUE)
          temptab<-table(temppar)
          if(length(temptab)==k&min(temptab)>=mingr&max(temptab)<=maxgr)find.new.par<-FALSE
        }
      }else{
        temppar<-NULL
        for(imode in 1:nmode){
          find.new.par<-TRUE
          while(find.new.par){
            if(ver!=4){
		itemppar<-integer(n[imode])

		if(ver==1){
			itemppar<-1:n[imode]%%k[imode]+1	
		}

		if(ver==2){
			itemppar[1:k[imode]]<-1:k[imode]
			itemppar[(k[imode]+1):n[imode]]<-k[imode]
		}

		if(ver==3){
			itemppar[1:k[imode]]<-1:k[imode]
			itemppar[(k[imode]+1):n[imode]]<-1+trunc(k[imode]*runif(n[imode]-k[imode]))
		}

		for(ii in n[imode]:2){
			jj<-trunc(ii*runif(1))
			itemppar[c(ii,jj)]<-itemppar[c(jj,ii)]
		}
            }else itemppar<-sample(1:k[imode],n[imode],replace=TRUE)
            temptab<-table(itemppar)
            if((length(temptab)==k[imode])&(min(temptab)>=mingr)&(max(temptab)<=maxgr))find.new.par<-FALSE
          }
          temppar<-c(temppar,list(itemppar))
        }
      }
      return(temppar)
}

