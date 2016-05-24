matrixBC <-
function(x,k,r,lambda,alpha,beta,nstart=20,Cs.init=NULL,Ds.init=NULL,max.iter=50,threshold=1e-4,Sigma.init=NULL,Delta.init=NULL,center=TRUE){
	
  if(is.null(Cs.init)){
	Cs<-kmeans(x,k,nstart=20)$cluster
  }
  else{
	Cs<-Cs.init
  }	
  if(is.null(Ds.init)){
    Ds <- kmeans(t(x), r,nstart=20)$cluster
  }
  else{
    Ds <- Ds.init
  }
  
  if(center==TRUE){
  	    mustemp <- mean(x)
    	x <- x-mustemp
  }
  
  cl <- match.call()
  Cslist<-list()
  Dslist<-list()

  if(is.null(Sigma.init)){
    Sigma<-diag(1,nrow=nrow(x),ncol=nrow(x))
    Delta<-diag(1,nrow=ncol(x),ncol=ncol(x))
    mus<-updateMusMatrix(x,Cs,Ds,lambda,Sigma,Delta)
    objs<-1e15
    improvement<-1e10
    i<-1 
    while(improvement>(threshold) && i<=max.iter){
      cat(i)
      
      mus<-updateMusMatrix(x,Cs,Ds,lambda,Sigma,Delta)

      Sigma<-updateSigma(Delta,mus[Cs,Ds],x,alpha,Cs,Ds)

      Delta<-updateDelta(Sigma,mus[Cs,Ds],x,beta,Cs,Ds)  

      Cs<-UpdateRowCluster(x,Sigma,Delta,mus[Cs,Ds],Cs,Ds)
      Cs<-ReNumberMatrix(Cs)
      Cslist[[i]]<-Cs

      mus<-updateMusMatrix(x,Cs,Ds,lambda,Sigma,Delta)

      Sigma<-updateSigma(Delta,mus[Cs,Ds],x,alpha,Cs,Ds)

      Delta<-updateDelta(Sigma,mus[Cs,Ds],x,beta,Cs,Ds)  

      Ds<-UpdateColumnCluster(x,Sigma,Delta,mus[Cs,Ds],Cs,Ds)
      Ds<-ReNumberMatrix(Ds)
      objs<-c(objs,MatrixObjective(x,mus[Cs,Ds],Cs,Ds,Sigma,Delta,lambda,alpha,beta))

      Dslist[[i]]<-Ds
      i<-i+1
      improvement<-abs(objs[i]-objs[i-1])
    }
  }  		  
  
############ If Sigma and Delta is initialized, we do not want to update them.    
  else{
	Sigma<-Sigma.init
    Delta<-Delta.init  	
 
    mus<-updateMusMatrix(x,Cs,Ds,lambda,Sigma,Delta)
    objs<-1e15
    improvement<-1e10
    i<-1
    
    while(improvement>threshold && i<=max.iter){
    mus<-updateMusMatrix(x,Cs,Ds,lambda,Sigma,Delta)
    Cs<-UpdateRowCluster(x,Sigma,Delta,mus[Cs,Ds],Cs,Ds)
    Cs<-ReNumberMatrix(Cs)
    Cslist[[i]]<-Cs

    mus<-updateMusMatrix(x,Cs,Ds,lambda,Sigma,Delta)

    Ds<-UpdateColumnCluster(x,Sigma,Delta,mus[Cs,Ds],Cs,Ds)
    Ds<-ReNumberMatrix(Ds)
    Dslist[[i]]<-Ds

    objs<-c(objs,MatrixObjective(x,mus[Cs,Ds],Cs,Ds,Sigma,Delta,lambda,alpha=0,beta=0))

    i<-i+1
    improvement<-abs(objs[i]-objs[i-1])
    }
  }	
  if(min(diff(objs))>0) print("Warning: objective values decreases")	
 
  if(center==TRUE){
		mus <- mus+mustemp
  }
 # if(improvement>1) stop("This is bad!!!")
 
  out <- list()
  class(out) <- "matrixBC"
  out$Cs <- Cslist[[i-1]]
  out$Ds <- Dslist[[i-1]]
  out$objs <- objs
  out$mus <- mus[out$Cs,out$Ds]
  out$Mus <- mus
  out$Sigma <- Sigma
  out$Delta <- Delta
  out$iteration <- i
  out$cl <- cl
    
  return(out)  
}
