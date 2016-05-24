nichepref <-
function (P, D = NULL, q = NULL, mode="multiple", Np = NULL, Nq = NULL, nboot = 1000, alpha=0.05) {
    if (!inherits(P, "data.frame")) stop("Non convenient dataframe for species resource use")
    if (!is.null(D)) {
        if (!inherits(D, "dist")) stop("Object of class 'dist' expected for distance")
        D <- as.matrix(D)
        if (ncol(P) != nrow(D)) stop("The number of columns in P must be equal to the number of items in D")
        D <- as.dist(D)
    }
    if(mode=="multiple" && !is.null(Np)) {
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
    	
    if(mode=="single" || !is.null(Np)) nc = 3
    else nc = 1

	  if(mode=="multiple") {
		  F=as.data.frame(matrix(NA, nrow=nrow(P),ncol= ncol(P)))
			names(F) = names(P)
 		  row.names(F)=row.names(P)
		  if(!is.null(Np))	 {
			 LF=as.data.frame(matrix(NA, nrow=nrow(P),ncol= ncol(P)))
			 row.names(LF)=row.names(P)
			 names(LF) = names(P)
			 UF=as.data.frame(matrix(NA, nrow=nrow(P),ncol= ncol(P)))
			 row.names(UF)=row.names(P)
			 names(UF) = names(P)
		  } 
	   for (i in 1:nrow(P)) {
	      pi = as.numeric(P[i,])
	    	 if(sum(is.na(getF(pi)))==0) {
	    		  	F[i,] = getF(pi,q)
	    		  	if(!is.null(Np)) {
	    		  	   #Generate bootstrap samples from multinomial distribution
	    		    	 BF = matrix(0,nrow=nboot, ncol=(ncol(P)))
	    		  	   bsamp = rmultinom(nboot,Np[i],getF(pi))
	    		  	   if(!is.null(Nq)) qsamp = rmultinom(nboot,Nq,q)
	    		  	   for(b in 1:nboot) {
	    		  	   	if(!is.null(Nq)) BF[b,] = getF(bsamp[,b],qsamp[,b])
	    		  	   	else BF[b,] = getF(bsamp[,b],q)
	    		  	   }
	    		  	   #Some NA may appear because of zeroes in qsamp
	    		  	   BF = BF[!is.na(rowSums(BF)),]
				    	 #Compute Bias-corrected percentile method (Manly 2007: pp52-56) for each dimension
	    		  	   for(j in 1:ncol(P)) {
		    		  	   z0 = qnorm(sum(BF[,j]<F[i,j])/nrow(BF))
	  	  		  	   	lj = floor(nrow(BF)*pnorm(2*z0+qnorm(alpha/2)))
	    		  	   	uj = floor(nrow(BF)*pnorm(2*z0+qnorm(1-(alpha/2))))
	    		  	   	if(lj > 0 && uj > 0) {
		 	   	  	  	 	sbf = sort(BF[,j])
		 	   	  	   	LF[i,j] = sbf[lj]
			   	  	   	UF[i,j] = sbf[uj]
		    	  	   	}
	    		  	   }
	    	  	  }
	    	  	}
	    }
   		if(!is.null(Np)) return(list(F=F, LF = LF, UF = UF))
	   else return(F)
	}
	else if(mode=="single") {
		 F=as.data.frame(matrix(NA, nrow=3,ncol= (ncol(P))))
		 row.names(F)=c("F", "F LC", "F UC")
		 names(F) = names(P)
		 
		 F[1,] = getF(colSums(P),q)

 	   #Generate bootstrap samples from multinomial distribution
	   BF = matrix(0,nrow=nboot, ncol=ncol(P))
	   if(!is.null(Nq)) qsamp = rmultinom(nboot,Nq,q)
 	   for(b in 1:nboot) {
	 	   psamp = colSums(P[sample(1:nrow(P),replace=TRUE),])
	    	 if(!is.null(Nq)) BF[b,] = getF(psamp,qsamp[,b])
	    	 else BF[b,] = getF(psamp,q)
	   }
	   #Some NA may appear because of zeroes in qsamp
	   BF = BF[!is.na(rowSums(BF)),]
		 #Compute Bias-corrected percentile method (Manly 2007: pp52-56) for each dimension
	   for(j in 1:ncol(P)) {
		    z0 = qnorm(sum(BF[,j]<F[1,j])/nrow(BF))
	  	  	  lj = floor(nrow(BF)*pnorm(2*z0+qnorm(alpha/2)))
	    		uj = floor(nrow(BF)*pnorm(2*z0+qnorm(1-(alpha/2))))
	    		if(lj > 0 && uj > 0) {
		 	   	sbf = sort(BF[,j])
		 	   	F[2,j] = sbf[lj]
			   	F[3,j] = sbf[uj]
		    	}
	   }
   	 return(F)
	}
}

