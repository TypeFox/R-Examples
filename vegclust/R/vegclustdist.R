vegclustdist <-
function(x,mobileMemb, fixedMemb=NULL, method="NC", m=2,dnoise=NULL, eta = NULL, alpha=0.001, iter.max=100, nstart=1, seeds = NULL, verbose=FALSE) {

#One run of vegclustdist   
vegclustonedist <-
function(d,mobileMemb, fixedMemb=NULL, method="NC", m=2,dnoise=NULL, eta = NULL, alpha=0.001, iter.max=100) {
   METHODS <- c("KM", "FCM", "PCM","NC","HNC" ,"KMdd","NCdd", "HNCdd", "FCMdd", "PCMdd")
   method <- match.arg(method, METHODS)
   if(method=="KM"||method=="KMdd") {
   	m=1.0
   	dnoise=NULL
   	eta=NULL
   }
   else if(method=="FCM"||method=="FCMdd") {
   	dnoise=NULL
   	eta=NULL
   }
   else if(method=="NC"||method=="NCdd") {
   	if(is.null(dnoise)) stop("Must provide a value for dnoise")
   	eta = NULL
   }
   else if(method=="HNC"||method=="HNCdd") {
     if(is.null(dnoise)) stop("Must provide a value for dnoise")
     eta = NULL
     m=1.0
   }
   else if(method=="PCM"||method=="PCMdd") {
   	if(is.null(eta)) stop("Must provide a vector of values for eta")
   	dnoise = NULL
   }
   
	d = as.matrix(d)
  #Number of objects
	n = nrow(d)
	
	#Sets the starting memberships for mobile clusters
	if(is.data.frame(mobileMemb) || is.matrix(mobileMemb)) {
		if(nrow(mobileMemb)!=ncol(d)) {
			stop("The number of rows in mobileMemb must be the same as the number rows and columns of d")
		}		
		u = as.matrix(mobileMemb)
	}
	else if(is.vector(mobileMemb) && length(mobileMemb)==1 && is.numeric(mobileMemb)) {
		s = sample(n, mobileMemb)
		u = matrix(0,n,length(s))
		for(k in 1:length(s)) {
			u[s[k],k]=1
		}
	} else if(is.vector(mobileMemb) && is.numeric(mobileMemb)) {
		u = matrix(0,n,length(mobileMemb))
		for(k in 1:length(mobileMemb)) {
			u[mobileMemb[k],k]=1
		}		
	}
	else {
		stop("Provide a number, a vector of seeds, or membership matrix for mobile clusters")
	}	
	kMov = ncol(u)
	
	#Sets the fixed cluster memberships
	if(!is.null(fixedMemb)) {
		if(is.data.frame(fixedMemb)) {
			kFix = ncol(fixedMemb)
			u = cbind(u,fixedMemb)
		}
		else if(!is.matrix(fixedMemb)) {
			stop("Fixed clusters must be specified as a matrix or a data frame")
		}	
		else {
			kFix = ncol(fixedMemb)
			u = cbind(u,fixedMemb)
		}
	} else {
		kFix = 0
	}
  #Define vector of medoids
  med = rep(NA,ncol(u))
   
	#Check possibilistic parameters
   if((method=="PCM"||method=="PCMdd") && length(eta)!=(kMov+kFix)) stop("Vector of reference distances (eta) must have a length equal to the number of clusters")

   #Add extra (noise) column for NC-related methods
   if(method=="NC"||method=="NCdd"||method=="HNC"||method=="HNCdd") {
   	u = cbind(u, vector("numeric",length=n))
  	}
  	uPrev = matrix(0,nrow=n,ncol=ncol(u))
	
   #Initialize squared distances to fixed centroids
   if(method=="KM"||method=="PCM"||method=="NC"|| method=="HNC"||method=="FCM") {     
     sqdist2cent = matrix(0,nrow=n,ncol=(kMov+kFix))
     if(kFix>0) {
   	  for(k in (kMov+1):(kMov+kFix)) {
			  vargeom = sum((u[,k]^m) %*% (d^2) %*% (u[,k]^m))/(2*sum(u[,k]^m)^2)
	  		  for(i in 1:n) {
	  			  sqdist2cent[i,k] = (sum((u[,k]^m)*(d[i,]^2))/sum(u[,k]^m))-vargeom
	  			  if(sqdist2cent[i,k]<0) sqdist2cent[i,k]=0
	  		  }
   	  }   	
     }
   } else { #Initialize distances to fixed medoids
     dist2med = matrix(0,nrow=n,ncol=(kMov+kFix))
     if(kFix>0) {
       for(k in (kMov+1):(kMov+kFix)) {
         #Determine medoid
         med[k] = which.min((u[,k]^m)%*%d)
         #Update distances to medoids
         dist2med[,k] = d[,med[k]] 
       }
     }
   }

	continue = TRUE
	iter = 1
   #iterates until no change in memberships
   while(continue) {
     #1. Update squared distances to centers (centroids for Euclidean-based methods and medoids for the others)
     if(method=="KM"||method=="PCM"||method=="NC"|| method=="HNC"||method=="FCM") {     
   	    vargeom = vector("numeric", kMov)
   	    for(k in 1:(kMov)) {
			    vargeom[k] = sum((u[,k]^m) %*% (d^2) %*% (u[,k]^m))/(2*sum(u[,k]^m)^2)
			    for(i in 1:n) {
			      sqdist2cent[i,k] = (sum((u[,k]^m)*(d[i,]^2))/sum(u[,k]^m))-vargeom[k]
			      if(sqdist2cent[i,k]<0) sqdist2cent[i,k]=0
			    }
   	    }
     } else{ 
       for(k in 1:kMov) {
         #Determine medoid
         med[k] = which.min((u[,k]^m)%*%d)
         dist2med[,k] = d[,med[k]]
       }
     }
     
     #2. compute membership to centroids for mobile and fixed clusters
     if (method=="KM") {
       minC<-apply(sqdist2cent,1,which.min)
       u[,] = 0
       for(k in 1:length(minC)) u[k,minC[k]] = 1.0
     } else if (method=="KMdd") {
         minC<-apply(dist2med,1,which.min)
         u[,] = 0
         for(k in 1:length(minC)) u[k,minC[k]] = 1.0
     } else if(method=="NC") {
       d2cm2<-cbind(sqdist2cent,dnoise^2)
       for(k in 1:ncol(d2cm2)) {
         a<-sweep(d2cm2,1,d2cm2[,k],"/")
         u[,k] = 1/rowSums(a^(-1/(m-1)))
       }
       u[d2cm2==0]=1
     } else if(method=="HNC") {
       d2cm<-cbind(sqdist2cent,dnoise^2)
       u[,] = 0
       minC<-apply(d2cm,1,which.min)
       for(k in 1:length(minC)) {
         u[k,minC[k]] = 1.0
       }
     } else if(method=="HNCdd") {
       d2cm<-cbind(dist2med,dnoise)
       u[,] = 0
       minC<-apply(d2cm,1,which.min)
       for(k in 1:length(minC)) {
         u[k,minC[k]] = 1.0
       }
     } else if(method=="NCdd") {
       d2cm<-cbind(dist2med,dnoise)
       for(k in 1:ncol(d2cm)) {
         a<-sweep(d2cm,1,d2cm[,k],"/")
         u[,k] = 1/rowSums(a^(-1/(m-1)))
       }
       u[d2cm==0]=1
     } else if (method=="FCM") {
       for(k in 1:ncol(sqdist2cent)) {
         a<-sweep(sqdist2cent,1,sqdist2cent[,k],"/")
         u[,k] = 1/rowSums(a^(-1/(m-1)))
       }
       u[sqdist2cent==0]=1
     } else if (method=="FCMdd") {
       d2cm<-dist2med
       for(k in 1:ncol(d2cm)) {
         a<-sweep(d2cm,1,d2cm[,k],"/")
         u[,k] = 1/rowSums(a^(-1/(m-1)))
       }
       u[dist2med ==0]=1
     } else if (method=="PCM") {
       for(k in 1:ncol(sqdist2cent)) u[,k] = 1/(1+((sqdist2cent[,k])/eta[k])^(1/(m-1)))
       u[dist2cent==0]=1
     } else if (method=="PCMdd") {
       for(k in 1:ncol(dist2med)) u[,k] = 1/(1+((dist2med[,k])/eta[k])^(1/(m-1)))
       u[dist2cent==0]=1
     } 	 	
       			
   	#Check for stopping
   	if(iter>2) {
   		continue = (max(abs(u-uPrev))>alpha) && (iter<=iter.max) && (max(abs(u-uPrev2))>alpha)
   	}   	
   	if(continue) {
	   	uPrev2 = uPrev
	   	uPrev = u	   	
	   	iter=iter+1
	   	if(verbose) cat(".")
   	}
   }
   if(method=="FCM" ||  method=="KM") functional = sum(sqdist2cent*(u^m))
   else if(method=="NC"||method=="HNC") functional = sum(sqdist2cent*(u[,-(kMov+kFix+1)]^m))+sum(dnoise^2*u[,kMov+kFix+1]^m)
   else if(method=="FCMdd"||  method=="KMdd") functional = sum(dist2med*(u^m))
   else if(method=="NCdd"||method=="HNCdd") functional = sum(dist2med*(u[,-(kMov+kFix+1)]^m))+sum(dnoise*u[,kMov+kFix+1]^m)
   else if(method=="PCM") {
    functional = 0
    for(k in 1:(kMov+kFix)) functional = functional+sum(sqdist2cent[,k]*(u[,k]^m))+sum(eta[k]*(1-u[,k])^m)
  } else if(method=="PCMdd") {
    functional = 0
    for(k in 1:(kMov+kFix)) functional = functional+sum(dist2med[,k]*(u[,k]^m))+sum(eta[k]*(1-u[,k])^m)
  } 
   if(verbose) cat(paste("\nIterations:", iter,"Functional: ", functional,"\n"))
   
   #Prepare output
  u = as.data.frame(u)   
  if(method=="KM"||method=="PCM"||method=="NC"|| method=="HNC"||method=="FCM") {     
    dist2cent = as.data.frame(sqrt(sqdist2cent))   
  } else {
    dist2cent = as.data.frame(dist2med)   
  }
	for(k in 1:kMov) {
		names(u)[k] = paste("M",k,sep="")
		names(dist2cent)[k] = paste("M",k,sep="")
	}
	if(kFix>1) {
		for(k in (kMov+1):(kMov+kFix)) {
			names(u)[k] = paste("F",k,sep="")
			names(dist2cent)[k] = paste("F",k,sep="")
		}
	}
	if(method=="NC"||method=="NCdd"||method=="HNC"||method=="HNCdd") names(u)[kMov+kFix+1] = "N"
	rownames(u) = rownames(d)
	rownames(dist2cent) = rownames(d)
	size = colSums(u[,1:(kMov+kFix)])
  if(method=="NC"||method=="FCM"||method=="KM"||method=="PCM"||method=="HNC") withinss = colSums((dist2cent^2)*(u[,1:(kMov+kFix)]^m))
  else withinss = colSums((dist2cent)*(u[,1:(kMov+kFix)]^m))
  #Return medoid indices as mobile or fixed centers
  mobileCenters = NULL
  fixedCenters = NULL
  if(method=="KMdd"||method=="FCMdd"||method=="NCdd"||method=="HNCdd"||method=="PCMdd") {
    mobileCenters = med[1:kMov]    
    if(kFix>0) fixedCenters = med[(kMov+1):(kMov+kFix)]
  }
  res = list(mode="dist", method=method, m = m, dnoise = dnoise,
             eta = eta, memb=u,
             mobileCenters=mobileCenters, fixedCenters= fixedCenters, 
             dist2clusters=dist2cent, withinss = withinss, 
             size=size, functional=functional)
  class(res)<-"vegclust"
	return(res)
}

	if(is.null(seeds)) seeds = 1:nrow(as.matrix(x))
   #If mobileCenters is a number and nstart>1 perform different random starts
	if(is.vector(mobileMemb) && length(mobileMemb)==1 && is.numeric(mobileMemb)) {
	   bestRun = vegclustonedist(x, mobileMemb=sample(seeds, mobileMemb), fixedMemb, method, m,dnoise, eta, alpha, iter.max)
		if(nstart>1) {
			for(i in 2:nstart) {
				run = vegclustonedist(x,mobileMemb=sample(seeds,mobileMemb), fixedMemb, method, m,dnoise,eta, alpha, iter.max)
				if(run$functional<bestRun$functional) {
					bestRun = run
				}
			}
		}		
		return(bestRun)
	} else { #Perform a single run
		return(vegclustonedist(x,mobileMemb, fixedMemb, method, m,dnoise, eta, alpha, iter.max))
	}
}

