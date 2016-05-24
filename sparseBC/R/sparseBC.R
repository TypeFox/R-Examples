sparseBC <-
function(x,k,r,lambda,nstart=20, Cs.init=NULL, Ds.init=NULL, max.iter=1000,threshold=1e-10,center=TRUE){
	
	
    if(is.null(Cs.init)){
      Cs <- (kmeans(x, k,nstart=nstart)$cluster)
    } else {
      Cs <- Cs.init
    }
    if(is.null(Ds.init)){
      Ds <- (kmeans(t(x), r,nstart=nstart)$cluster)
    } else {
      Ds <- Ds.init
    }
    	
    if(center==TRUE)
    {
    	mustemp <- mean(x)
    	x <- x-mustemp
    }
    
    cl <- match.call()
	mus <- UpdateMus(x, Cs, Ds,lambda=lambda)
	objs <- 1e15
	improvement <- 1e10
	i<-1
	
	while(improvement>(threshold) && i<=max.iter){
		Cs <- UpdateClusters(x,mus,Cs,Ds)
		objs <- c(objs, Objective(x, mus,Cs,Ds,lambda=lambda))
		Cs<-ReNumber(Cs)
		
		mus <- UpdateMus(x, Cs, Ds, lambda=lambda)
		objs <- c(objs, Objective(x, mus,Cs,Ds,lambda=lambda))
		
		Ds <- UpdateClusters(t(x),t(mus),Ds,Cs)
		objs <- c(objs, Objective(x, mus,Cs,Ds,lambda=lambda))
		Ds<-ReNumber(Ds)
		
		mus <- UpdateMus(x, Cs, Ds,lambda=lambda)
		objs <- c(objs, Objective(x, mus,Cs,Ds,lambda=lambda))
		
		improvement <- abs(objs[length(objs)]-objs[length(objs)-4])/abs(objs[length(objs)-4])	
		i<-i+1
	}
	
	if(i > max.iter){
		warning("The algorithm has not converged by the specified maximum number of iteration")
	}
	
	if(center==TRUE){
		mus <- mus+mustemp
	}
	
 # return output
 	out <- list()
 	class(out) <- "sparseBC"	
 	out$Cs <- Cs
 	out$Ds <- Ds
 	out$objs <- objs
 	out$mus <- mus[Cs,Ds]
 	out$Mus <- mus
 	out$iteration <- i 	 	 		
	out$cl <- cl
	
	return(out)
}








