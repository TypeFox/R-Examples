nichecentroid <-
function (P, D = NULL, q = NULL, mode="multiple", Np = NULL, Nq = NULL, nboot = 1000, alpha=0.05) {
    if (!inherits(P, "data.frame")) stop("Non convenient dataframe for species resource use")
    if (!is.null(D)) {
        if (!inherits(D, "dist")) stop("Object of class 'dist' expected for distance")
        D <- as.matrix(D)
        if (ncol(P) != nrow(D)) stop("The number of columns in P must be equal to the number of items in D")
        D <- as.dist(D)
    }
    if(!is.null(Np)) {
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
	 
    	# Returns preference from a resource use vector (considering resource availability in desired)
    	getF<-function(p,q=NULL) {
    		if(!is.null(q)) {
    			a = p/q
    			return(a/sum(a))
    		} else { #Just to check that we have proportions
    			return(p/sum(p))
    		}
    	}
    	
    		 
	 #Computes metric MDS
	 cmd = cmdscale(D,eig=TRUE,k= ncol(P)-1)
    
   if(mode=="multiple") {
		 C=as.data.frame(matrix(NA, nrow=nrow(P),ncol= ncol(P)-1))
		 row.names(C)=row.names(P)
		 if(!is.null(Np))	 {
			 LC=as.data.frame(matrix(NA, nrow=nrow(P),ncol= ncol(P)-1))
			 row.names(LC)=row.names(P)
			 UC=as.data.frame(matrix(NA, nrow=nrow(P),ncol= ncol(P)-1))
			 row.names(UC)=row.names(P)
		 } 
	   for (i in 1:nrow(P)) {
	      pi = as.numeric(P[i,])
	    	 if(sum(is.na(getF(pi)))==0) {
	    		  	C[i,] = (getF(pi,q)%*%cmd$points)/sum(getF(pi,q))
					if(!is.null(Np)) {
 	   		  	   #Generate bootstrap samples from multinomial distribution
	    		    	 BC = matrix(0,nrow=nboot, ncol=(ncol(P)-1))
 	   		  	   bsamp = rmultinom(nboot,Np[i],getF(pi))
	    		  	   if(!is.null(Nq)) qsamp = rmultinom(nboot,Nq,q)
 	   		  	   for(b in 1:nboot) {
	    		  	   	if(!is.null(Nq)) BC[b,] = (getF(bsamp[,b],qsamp[,b])%*%cmd$points)/sum(getF(bsamp[,b],qsamp[,b]))
	    		  	   	else BC[b,] = (getF(bsamp[,b],q)%*%cmd$points)/sum(getF(bsamp[,b],q))
	    		  	   }
	    		  	   #Some NA may appear because of zeroes in qsamp
	    		  	   BC = BC[!is.na(rowSums(BC)),]
				    	 #Compute Bias-corrected percentile method (Manly 2007: pp52-56) for each dimension
 	   		  	   for(j in 1:(ncol(P)-1)) {
		    		  	   z0 = qnorm(sum(BC[,j]<C[i,j])/nrow(BC))
  		  		  	   	lj = floor(nrow(BC)*pnorm(2*z0+qnorm(alpha/2)))
 	   		  	   	uj = floor(nrow(BC)*pnorm(2*z0+qnorm(1-(alpha/2))))
 	   		  	   	if(lj > 0 && uj > 0) {
		 	   	  	  	 	sbc = sort(BC[,j])
		 	   	  	   	LC[i,j] = sbc[lj]
			   	  	   	UC[i,j] = sbc[uj]
		    	  	   	}
 	   		  	   }
 	   	  	  }
 	   	  	}
 	   }
	   if(!is.null(Np)) return(list(C=C, LC = LC, UC = UC))
  		 else return(C)
   }
   else if(mode=="single") {
		 C=as.data.frame(matrix(NA, nrow=3,ncol= ncol(P)-1))
		 row.names(C)=c("Centroid", "Centroid LC", "Centroid UC")
	   C[1,] = (getF(colSums(P),q)%*%cmd$points)/sum(getF(colSums(P),q))
 	   #Generate bootstrap samples from multinomial distribution
	   BC = matrix(0,nrow=nboot, ncol=(ncol(P)-1))
	   if(!is.null(Nq)) qsamp = rmultinom(nboot,Nq,q)
 	   for(b in 1:nboot) {
	 	   psamp = colSums(P[sample(1:nrow(P),replace=TRUE),])
	    	 if(!is.null(Nq)) BC[b,] = (getF(psamp,qsamp[,b])%*%cmd$points)/sum(getF(psamp,qsamp[,b]))
	    	 else BC[b,] = (getF(psamp,q)%*%cmd$points)/sum(getF(psamp,q))
	    }
	    #Some NA may appear because of zeroes in qsamp
	    BC = BC[!is.na(rowSums(BC)),]
		  #Compute Bias-corrected percentile method (Manly 2007: pp52-56) for each dimension
 	    for(j in 1:(ncol(P)-1)) {
		     z0 = qnorm(sum(BC[,j]<C[1,j])/nrow(BC))
  		  		lj = floor(nrow(BC)*pnorm(2*z0+qnorm(alpha/2)))
 	   		uj = floor(nrow(BC)*pnorm(2*z0+qnorm(1-(alpha/2))))
 	   		if(lj > 0 && uj > 0) {
		 	   	 sbc = sort(BC[,j])
		 	   	 C[2,j] = sbc[lj]
			   	 C[3,j] = sbc[uj]
		    	 }
 	   	}
   	  return(C)
   } 
    
}

