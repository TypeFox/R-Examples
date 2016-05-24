
####################################################################################################################
## Author: Gro Nilsen
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the clusterGenomics package
## Reference: "Identifying clusters in genomics data by recursive partitioning", Nilsen et al. (2013, preprint)
####################################################################################################################


### Requires: ##
## clusteringFunctions.r
	
### Required by: ##
## part.r

#gap finds the optimal number of clusters based on the "gap statistic" (Tibshirani et al. 2001) 

gap <- function(X,Kmax=10,B=100,ref.gen="PC",cl.lab=NULL,...){
	
  #Default ... values:
  default.par <- list(dist.method="euclidean",cl.method="hclust",linkage="average",cor.method="pearson",nstart=10)
  #Check for user modifications
  #Note: call could come from part in which case ... is a list called fixed.par, or it could be a independent call in which case ... could contain several parameters which must be converted to list
  if(hasArg(fixed.par)){
    fixed.par <- modifyList(default.par,as.list(list(...)$fixed.par))
  }else{
    fixed.par <- modifyList(default.par,list(...))    
  }
  
  n <- nrow(X)
	
  if(n<=Kmax){
    if(fixed.par$cl.method=="hclust"){
      Kmax = n
    }else{
      Kmax = n-1   #kmeans does not work if K=n
    }
  }
  	
	#Calculate distances:
	dX <- getDist(X,dist.method=fixed.par$dist.method,cor.method=fixed.par$cor.method)
	
	if(is.null(cl.lab)){
	 cl.lab <- findPartition(X=X,Kmax=Kmax,dX=dX,fixed.par)
	}
	
	#Find W: vector containing W_K for all choices of K
	W <- findW(dX=dX,K=Kmax,cl.lab=cl.lab)

		
	#Generate reference data sets and find Wb:
	Wb <- getReferenceW(X,Kmax,B,ref.gen,fixed.par)

	
	#Calculate gap statistic:
	L <- apply(log(Wb),1,mean)
  gap <- L - log(W)

	sdk <- apply(log(Wb),1,sd)*sqrt((B-1)/B)  #multiply by the last expression to get 1/B (as in original article) instead of 1/(B-1) (normal calculation)
	sk <- sqrt(1+(1/B))*sdk

	#Calculate gap-criterion for each k:
	diff <- gap[1:Kmax-1] - (gap[2:Kmax] - sk[2:Kmax])

	#Select the first k where the criterion is positive:
	kvec <- 1:Kmax
	posDiff <- diff[diff>=0]
	if(length(posDiff)==0){
		hatK <- 1  #no k satisfies criterion, define 1 cluster e.g. to make part run..
	}else{ 
		hatK <- kvec[diff>=0][1]
	}

	if(all(is.na(diff))){
		hatK <- 1   #when used with part; can happen if sub-cluster consists of only 1 or 2 objects
	}
	
  #Get labels for best partition:
	lab.hatK <- cl.lab[[hatK]]
		
	return(list(hatK=hatK,lab.hatK=lab.hatK,gap=gap,sk=sk,W=W))
}



## helper-functions only used by GAP:

#Simulate from a uniform distribution according to min and max of the feature vector:
sim <- function(Xcol) {
	min <- min(Xcol)
	max <- max(Xcol)
	U <- runif(length(Xcol),min,max)
	return(U)
}


#Calulates the total within-cluster dispersion (W_K) for a k=1,...,K:
findW <- function(dX,K,cl.lab){

	W <- rep(0,K)
	#n <- nrow(X)
	n <- nrow(as.matrix(dX))
	
	#Calculate W for k=1:
	W[1] <- sum(as.matrix(dX))/(2*n)

	k <- 2
	while(k<=K){
	  labX <- cl.lab[[k]]
		
		for(i in 1:k){
			#Sum of within-cluster dispersion:
			d.k <- as.matrix(dX)[labX==i,labX==i]
			D.k <- sum(d.k)
				
			nk <- nrow(as.matrix(d.k))
			W[k] <- W[k] + D.k/(2*nk)
			
		}#endfor
		k <- k+1
	}#endwhile
	return(W)

}#endfunction



getReferenceW <- function(X,Kmax,B,ref.gen,...)	{
  arg <- as.list(...)
	#Generate reference data sets and find Wb:
	Wb <- matrix(0,nrow=Kmax,ncol=B)
	if(ref.gen=="PC"){
    #Transform data using svd:
    m <- apply(X,2,mean,na.rm=T)   #First columncenter X:
    Xc <- sweep(X,2,m)
    #SVD:
    s <- svd(Xc)	
	  newX <- Xc%*%s$v
  }
	for(b in 1:B){
		if(ref.gen=="PC"){
      U <- apply(newX,2,sim)
		  Z1 <- U%*%t(s$v)      #backtransform
		  #Add mean
		  Z <- sweep(Z1,2,m,FUN="+")
		}else{
		  Z <- apply(X,2,sim)
		}
		#Calculate distances
		dZ <- getDist(Z,dist.method=arg$dist.method,cor.method=arg$cor.method)
		#Cluster reference data 
		clW.lab <- findPartition(X=Z,Kmax=Kmax,dX=dZ,arg)
		#Calculate Wb_K for all values of K
    Wb[,b] <- findW(dX=dZ,K=Kmax,cl.lab=clW.lab)
	}#endfor
	
	return(Wb)
	
}	