nichevar <-
function (P, D = NULL, q = NULL, mode="multiple", Np = NULL, Nq = NULL, nboot = 1000, alpha=0.05) {
    if (!inherits(P, "data.frame")) stop("Non convenient dataframe for species resource use")
    if (!is.null(D)) {
        if (!inherits(D, "dist")) stop("Object of class 'dist' expected for distance")
        D <- as.matrix(D)
        if (ncol(P) != nrow(D)) stop("The number of columns in P must be equal to the number of items in D")
        D <- as.dist(D)
    }
    if(!is.null(Np) && mode=="multiple") {
    	   if(length(Np)!=nrow(P)) stop("The number of items in Np must be equal to the number of rows in P")
    }
    if(!is.null(q)) {
    	   if(length(q)!=ncol(P)) stop("The number of items in q must be equal to the number of columns in P")
    	   q = q/sum(q) #Just to check that we have proportions
    } else {
        #If no availability is provided, then all resources are equally available
    	   q = rep(1/ncol(P),ncol(P))
    }
    	
    #If no distance matrix is provided, the distance between resources is assumed to be maximum
    if (is.null(D)) D <- as.dist((matrix(1, ncol(P), ncol(P)) - diag(rep(1, ncol(P)))))
    
    # Computes the niche breadth from the resource preferences of the target and the resource relationships
    nichevar1<-function(f, D) {
       if (is.na(sum(f))) v <- NA
       else if (sum(f) < 1e-16) v <- 0
       else v <- (f %*% (as.matrix(D)^2) %*% f)/(2*(sum(f)^2))
   		return(v)
    	}
    	
    	# Returns preference from a resource use vector (considering resource availability in desired)
    	getF<-function(p,q=NULL) {
    		if(!is.null(q)) {
    			a = p/q
    			return(a/sum(a))
    		} else { #Just to check that we have proportions
    			return(p/sum(p))
    		}
    	}
    	
    if(!is.null(Np) || mode=="single") nc = 3
    else nc = 1
    
    #Rows in P are different niches
    if(mode=="multiple") {
	    B <- as.data.frame(matrix(0,nrow=nrow(P), ncol=nc))
 	   rownames(B) <- row.names(P)
 	   for (i in 1:nrow(P)) {
 	   	  pi = as.numeric(P[i,])
 	   	  B[i,1] = nichevar1(getF(pi,q), D)
 	   	  
 	   	  if(!is.null(Np)) {
  	  		  		BB = vector("numeric",length=nboot)
 	   		  	if(sum(is.na(getF(pi)))==0) {
	    		  	   #Generate bootstrap samples from multinomial distribution
 	   		  	   psamp = rmultinom(nboot,Np[i],getF(pi))
 	   		  	   if(!is.null(Nq)) qsamp = rmultinom(nboot,Nq,q)
 	   		  	   for(b in 1:nboot) {
 	   		  	   	 if(!is.null(Nq)) BB[b] = nichevar1(getF(psamp[,b],qsamp[,b]),D)
 	   		  	   	 else BB[b] = nichevar1(getF(psamp[,b],q),D)
 	   		  	   }
  	  		  	   	#Some NA may appear because of zeroes in qsamp
	    		  	   BB = BB[!is.na(BB)]
				    	 #Compute Bias-corrected percentile method (Manly 2007: pp52-56)
 	   		  	   z0 = qnorm(sum(BB<B[i,1])/length(BB))
 	   		  	   lj = floor(length(BB)*pnorm(2*z0+qnorm(alpha/2)))
 	   		  	   uj = floor(length(BB)*pnorm(2*z0+qnorm(1-(alpha/2))))
 	   		  	   if(lj > 0 && uj > 0) {
		 	   	  	   sbb = sort(BB)
		 	   	  	   B[i,2] = sbb[lj]
			   	  	   B[i,3] = sbb[uj]
		    	  	   }
 	   	  	  }
 	   	  	}
 	   }
 	  }
    #Rows in P are observations
    else if(mode=="single") {
	   B <- as.data.frame(matrix(0,nrow=1, ncol=nc))
 	   rownames(B) <- "Niche"
 	   B[1,1] = nichevar1(getF(colSums(P),q), D)
 	   	  
  	  	 BB = vector("numeric",length=nboot)
 	   if(!is.null(Nq)) qsamp = rmultinom(nboot,Nq,q)
  	  	 for(b in 1:nboot) {
	 	   psamp = colSums(P[sample(1:nrow(P),replace=TRUE),])
	 	   if(!is.null(Nq)) {
 		   	 BB[b] = nichevar1(getF(psamp,qsamp[b]),D)
 		   } else {
 		   	 BB[b] = nichevar1(getF(psamp,q),D)
 		   }
  	  	 }
  	  	#Some NA may appear because of zeroes in qsamp
	   BB = BB[!is.na(BB)]
	   #Compute Bias-corrected percentile method (Manly 2007: pp52-56)
 	   z0 = qnorm(sum(BB<B[1,1])/length(BB))
 	   lj = floor(length(BB)*pnorm(2*z0+qnorm(alpha/2)))
 	   uj = floor(length(BB)*pnorm(2*z0+qnorm(1-(alpha/2))))
 	   if(lj > 0 && uj > 0) {
		 	  sbb = sort(BB)
		 	  B[1,2] = sbb[lj]
			  B[1,3] = sbb[uj]
		 }
 	  } 
 	  if(nc==1) names(B) <- "B"
 	  else names(B) <- c("B","LC", "UC")
    return(B)
}

