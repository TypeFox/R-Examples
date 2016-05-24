sparseBC.choosekr <-
function(x,k,r,lambda,percent=0.1,trace=FALSE){
	if((1%%percent)!=0) stop("1 must be divisible by the specified percentage")
	if(percent<=0) stop("percentage cannot be less than or equal to 0")
	if(percent>=1) stop("percentage cannot be larger or equal to 1")
	if(sum(diff(k)<=0)>0 || sum(diff(r)<=0)>0) stop("k and r has to be an increasing sequence.  Please sort k and r before using the function")

# Initialize a couple of variables
  	miss<-sample(1:(nrow(x)*ncol(x)),nrow(x)*ncol(x),replace=FALSE)
  	numberoftimes<-1/percent        
  	allresults<-array(NA, dim=c(numberoftimes, length(k), length(r)))

          Cs.init <- matrix(NA, nrow=nrow(x), ncol=length(k))
          for(i in 1:length(k)){
            Cs.init[,i] <- kmeans(x, k[i], nstart=20)$cluster
          }
          Ds.init <- matrix(NA, nrow=ncol(x), ncol=length(r))
          for(j in 1:length(r)){
            Ds.init[,j] <- kmeans(t(x), r[j], nstart=20)$cluster
          }
        	
	for(i in 1:numberoftimes){
        if(trace==TRUE) cat("Iteration", i, fill=TRUE)
		xmiss<-x
        missing <- miss[1:(nrow(x)*ncol(x))/numberoftimes]
		xmiss[missing]<-NA
#		xmiss[missing]<- (outer(apply(xmiss, 1, mean, na.rm=TRUE), apply(xmiss, 2, mean, na.rm=TRUE), "+")-mean(xmiss, na.rm=TRUE))[missing]
        xmiss[missing] <- mean(xmiss, na.rm=TRUE)
		
# Perform biclustering for various k and r specified and calculate MSE          
		for(a in 1:length(k)){
			for(b in 1:length(r)){
				res<-sparseBC(xmiss,k[a],r[b],lambda=lambda, Cs.init=Cs.init[,a], Ds.init=Ds.init[,b])$mus
                allresults[i,a,b] <- sum((x[missing]-res[missing])^2)
			}
		}
		miss<-miss[-1:-(dim(x)[1]*dim(x)[2]/numberoftimes)]
	}
	
    results.se <- apply(allresults, c(2,3), sd)/sqrt(numberoftimes)
    results.mean <- apply(allresults, c(2,3), mean)
        
    IndicatorMatrix <- 1*(results.mean[1:(length(k)-1), 1:(length(r)-1)] <=  results.mean[2:length(k), 2:length(r)]+results.se[2:length(k), 2:length(r)])
    
    
    if(max(IndicatorMatrix)==0) return(list(bestK=max(k), bestR=max(r))) # It seems best to have the largest value of (k,r)...
   
    RowIndexPlusColIndex <- outer(k[-length(k)],r[-length(r)],"*")
    smallestIndicatorTrue <- min(RowIndexPlusColIndex[IndicatorMatrix==TRUE])
    out <- which(IndicatorMatrix==TRUE & RowIndexPlusColIndex==smallestIndicatorTrue, arr.ind=TRUE)
    
# I want to return the true number of cluster k and r
    
    out<-matrix(c(k[out[,1]],r[out[,2]]),nrow(out),ncol(out))
    temprow<-NULL
    tempcol<-NULL
    for(i in 1:length(k)){
    	temprow<-c(temprow,paste("K = ",k[i],sep=""))
    }
    for(i in 1:length(r)){
    	tempcol<-c(tempcol,paste("R = ",r[i],sep=""))
    }
        
    rownames(results.se)<-temprow
    colnames(results.se)<-tempcol
    
    
    rownames(results.mean)<-temprow
    colnames(results.mean)<-tempcol
    
  return(list(estimated_kr=out, results.se=results.se, results.mean=results.mean))

}
