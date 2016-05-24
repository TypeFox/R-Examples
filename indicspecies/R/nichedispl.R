nichedispl <-
function (P1,P2 = NULL, D = NULL, q1 = NULL, q2 = NULL, mode="multiple", Np1 = NULL, Np2=NULL, Nq1 = NULL, Nq2 = NULL,nboot = 1000, alpha = 0.05) {
  MODES <- c("single", "multiple", "pairwise")
  mode <- match.arg(mode, MODES)
  if(mode =="multiple" || mode =="single") {
    if(is.null(P2)) stop("P2 cannot be null in mode 'multiple' or 'single'")
    if (!inherits(P1, "data.frame") || !inherits(P2, "data.frame")) stop("P1 and P2 should be dataframes")
    if (mode == "multiple" && (nrow(P1) != nrow(P2))) stop("Resource use dataframes do not have the same number of rows")
    if (ncol(P1) != ncol(P2)) stop("Resource use dataframes do not have the same number of columns (resources)")  
    if (!is.null(Np1) && !is.null(Np2)) {
      if (!inherits(Np1, "numeric") || !inherits(Np2, "numeric")) stop("Np1 and Np2 should be numeric vectors")
      Np1 = as.integer(Np1)
      Np2 = as.integer(Np2)
    }
    if((!is.null(Nq1) && is.null(Nq2)) || (is.null(Nq1) && !is.null(Nq2))) stop("Nq1 and Nq2 should be both either NULL or contain numeric values")
    if (!is.null(Nq1) && !is.null(Nq2)) {
      if (!inherits(Nq1, "numeric") || !inherits(Nq2, "numeric")) stop("Nq1 and Nq2 should be numeric")
      Nq1 = as.integer(Nq1)
      Nq2 = as.integer(Nq2)
    }
    if (!is.null(D)) {
      if (!inherits(D, "dist")) 
        stop("Object of class 'dist' expected for distance")
      D <- as.matrix(D)
      if (ncol(P1) != nrow(D)) stop("The number of columns in P1 must be equal to the number of items in D")
      if (ncol(P2) != nrow(D)) stop("The number of columns in P2 must be equal to the number of items in D")
      D <- as.dist(D)
    }
  } else if(mode=="pairwise") {
    if (!inherits(P1, "data.frame")) stop("P1 should be a dataframe")
    if (!is.null(D)) {
      if (!inherits(D, "dist")) stop("Object of class 'dist' expected for distance")
      D <- as.matrix(D)
      if (ncol(P1) != nrow(D)) stop("The number of columns in P1 must be equal to the number of items in D")
      D <- as.dist(D)
    }    
  }
  #Check 'q1' the availability of resources during the first 'season'
  if (!is.null(q1)) {
    if (length(q1) != ncol(P1)) stop("The number of items in q must be equal to the number of columns in P1 and P2")
    q1 = q1/sum(q1)
  } else {
    q1 = rep(1/ncol(P1), ncol(P1))
  }
  #Check 'q2' the availability of resources during the second 'season'
  if (!is.null(q2)) {
    if (length(q2) != ncol(P2)) stop("The number of items in q must be equal to the number of columns in P1 and P2")
    q2 = q2/sum(q2)
  } else {
    q2 = rep(1/ncol(P2), ncol(P2))
  }
  #If no distance matrix was supplied, generate one (equidistant resources)
  if (is.null(D)) D <- as.dist((matrix(1, ncol(P1), ncol(P1)) - diag(rep(1, ncol(P1)))))
  
  # INTERNAL FUNCTIONS
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
  #Computes niche centroid distance from two distribution of resource preferences and the resource relationships
  disp1<-function(f1,f2, D) {
    if (is.na(sum(f1)) || is.na(sum(f2))) cd <- NA
    else if (sum(f1) < 1e-16 || sum(f2) < 1e-16) cd <- 0
    else {
      cd <- sqrt(((f1 %*% (as.matrix(D)^2) %*% f2)/(sum(f1)*sum(f2))) - nichevar1(f1,D)-nichevar1(f2,D))
    }
    return(cd)
  }
    	  
    
  #Calculate displacement between each row of P1 and the corresponding row in P2
  if(mode=="multiple") {
    if((!is.null(Np1) && !is.null(Np2)))  nc = 3
    else nc = 1
    CD <- as.data.frame(matrix(NA,nrow=nrow(P1), ncol=nc))
	  for(i in 1:nrow(P1)) rownames(CD)[i] <- paste(row.names(P1)[i],"vs",row.names(P2)[i])
	  for (i in 1:nrow(P1)) {
	    pi1 = as.numeric(P1[i,])
	    pi2 = as.numeric(P2[i,])
	    CD[i,1] = disp1(getF(pi1,q1),getF(pi2,q2),D) 
	    if(!is.null(Np1) && !is.null(Np2)) {
	    	BCD = vector("numeric",length=nboot)
	    	#Generate bootstrap samples from multinomial distribution
 	   	  if(!is.na(sum(pi1)) && !is.na(sum(pi2))) {
	 	   	  bsamp1 = rmultinom(nboot,Np1[i],getF(pi1))
	 	   	  bsamp2 = rmultinom(nboot,Np2[i],getF(pi2))
	    		if(!is.null(Nq1) && !is.null(Nq2)) {
            qsamp1 = rmultinom(nboot,Nq1,q1)
            qsamp2 = rmultinom(nboot,Nq2,q2)
	    		}
	    		for(b in 1:nboot) {
	 	   	    if(!is.null(Np1) && !is.null(Np2)) BCD[b] = disp1(getF(bsamp1[,b], qsamp1[,b]),getF(bsamp2[,b], qsamp2[,b]),D)
	 	   	  	else BCD[b] = disp1(getF(bsamp1[,b], q1),getF(bsamp2[,b], q2),D)
	 	   	  }
	    		#Some NA may appear because of zeroes in qsamp
          BCD = BCD[!is.na(BCD)]
				  #Compute Bias-corrected percentile method (Manly 2007: pp52-56)
	 	   	  z0 = qnorm(sum(BCD<CD[i,1])/length(BCD))
	 	   	  lj = floor(length(BCD)*pnorm(2*z0+qnorm(alpha/2)))
	 	   	  uj = floor(length(BCD)*pnorm(2*z0+qnorm(1-(alpha/2))))
	 	   	  if(lj > 0 && uj > 0 && lj!=uj) {
			    	sbcd = sort(BCD)
            CD[i,2] = sbcd[lj]
			   	  CD[i,3] = sbcd[uj]
		    	}
 	   	  }
 	   	}
 	  }
    if(nc==1) names(CD) <- "CD"
    else names(CD) <- c("CD","LC", "UC")
    return(CD)    
 	}
  
  #P1 are different resource use assessments of a single entity (the same for P2)
  if(mode=="single") {
	  CD <- as.data.frame(matrix(NA,nrow=1, ncol=3))
	  rownames(CD) <- "Centr.Displ."
	  CD[1, 1] <- disp1(getF(colSums(P1, na.rm=TRUE), q),getF(colSums(P2, na.rm=TRUE), q),D)
	  BCD = vector("numeric",length=nboot)
	  if(!is.null(Nq1) && !is.null(Nq2)) {
	    qsamp1 = rmultinom(nboot,Nq1,q1)
	    qsamp2 = rmultinom(nboot,Nq2,q2)
	  }
	  for (b in 1:nboot) {
		 	p1samp = colSums(P1[sample(1:nrow(P1),replace=TRUE),], na.rm=TRUE)
		 	p2samp = colSums(P2[sample(1:nrow(P2),replace=TRUE),], na.rm=TRUE)
 		  if(!is.null(Nq1) && !is.null(Nq2)) BCD[b] = disp1(getF(p1samp, qsamp1[,b]),getF(p2samp, qsamp2[,b]),D)
 		  else BCD[b] = disp1(getF(p1samp, q1),getF(p2samp, q2),D)
	 	}
 	  #Some NA may appear because of zeroes in qsamp
	  BCD = BCD[!is.na(BCD)]
		#Compute Bias-corrected percentile method (Manly 2007: pp52-56)
 		z0 = qnorm(sum(BCD<CD[1,1])/length(BCD))
 		lj = floor(length(BCD)*pnorm(2*z0+qnorm(alpha/2)))
 		uj = floor(length(BCD)*pnorm(2*z0+qnorm(1-(alpha/2))))
 		if(lj > 0 && uj > 0 && lj!=uj) {
		  sbcd = sort(BCD)
		  CD[1,2] = sbcd[lj]
		  CD[1,3] = sbcd[uj]
		}
	  names(CD) <- c("CD", "LC", "UC")
	  return(CD)
  }
  
  #Only P1 is used. Calculate overlap between pairs of rows in P1. No confidence intervals are calculated
  if(mode=="pairwise") {
    CD <- matrix(0, nrow = nrow(P1), ncol = nrow(P1))
    rownames(CD)<-rownames(P1)
    colnames(CD)<-rownames(P1)
    if(!is.null(Np1)) {
      LC <- CD
      UC <- CD
    }
    for (i in 1:(nrow(P1)-1)) {
      for (j in (i+1):nrow(P1)) {
        pi = as.numeric(P1[i, ])
        pj = as.numeric(P1[j, ])
        CD[i,j] <- disp1(getF(pi, q), getF(pj, q), D)
        CD[j,i] <- CD[i,j]
        if (!is.null(Np1)) {
          BCD = vector("numeric", length = nboot)
          if (!is.na(sum(pi)) && !is.na(sum(pj))) {
            bsampi = rmultinom(nboot, Np1[i], getF(pi))
            bsampj = rmultinom(nboot, Np1[j], getF(pj))
            if (!is.null(Nq1)) qsamp1 = rmultinom(nboot, Nq1, q1)
            for (b in 1:nboot) {
              if (!is.null(Nq1)) BCD[b] = disp1(getF(bsampi[, b], qsamp1[, b]), getF(bsampj[, b], qsamp1[, b]), D)
              else BCD[b] = disp1(getF(bsampi[, b], q), getF(bsampj[, b], q), D)
            }
            BCD = BCD[!is.na(BCD)]
            z0 = qnorm(sum(BCD < CD[i,j])/length(BCD))
            lj = floor(length(BCD) * pnorm(2 * z0 + qnorm(alpha/2)))
            uj = floor(length(BCD) * pnorm(2 * z0 + qnorm(1 - (alpha/2))))
            if (lj > 0 && uj > 0 && lj != uj) {
              sbcd = sort(BCD)
              LC[i, j] = LC[j, i] =sbcd[lj]
              UC[i, j] = UC[j, i] =sbcd[uj]
            }
          }
        }
      }
    }
    if(!is.null(Np1)) {
      return(list(CD=CD, LC=LC, UC = UC))
    }
    else return(CD)
  } 	  	
}

