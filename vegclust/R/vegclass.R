vegclass<-function(y, x) {
	if(!inherits(y,"vegclust")) stop("y must be a vegclust object")
	

	if(y$mode=="raw"){
	    if(!is.null(y$fixedCenters)) centers = as.data.frame(rbind(y$mobileCenters,y$fixedCenters))
		  else centers = as.data.frame(y$mobileCenters)
		  if(length(names(x))!=length(names(centers)) || sum(names(x)==names(centers))<ncol(x)) {
				  c = conformveg(x,centers)
				  x = as.matrix(c$x)
				  centers = as.matrix(c$y)
		  } else {
			    x = as.matrix(x)
			    centers = as.matrix(centers)
		  }
		  k = nrow(centers)
	} else { 
		memb = as.matrix(y$memb)
    x = as.matrix(x) #x contains the distance from new objects to old ones
		d2cl = as.matrix(y$dist2clusters)
		k = ncol(d2cl)
	}

	m = y$m
  dnoise = y$dnoise
  eta = y$eta
  method = y$method
	n = nrow(x)
		

  if(method=="NC"||method=="NCdd"||method=="HNC"||method=="HNCdd") {
   	u = matrix(0,nrow=n,ncol=(k+1))
  } else {
   	u = matrix(0,nrow=n,ncol=k)
  }
	

	dist2cent = matrix(0,nrow=nrow(x),ncol=k)

	#1. compute distance to centers (fixed and mobile)
	if(y$mode=="raw") {
	  	for(i in 1:k) {
		   dist2cent[,i] = sqrt(rowSums(sweep(x,2,centers[i,],"-")^2))
		}
	} else { #distance mode
	  if(method=="KM"||method=="PCM"||method=="NC"|| method=="HNC"||method=="FCM") {     
  		for(i in 1:k) {
	  		vargeom = sum((memb[,i]^m) %*% (d2cl[,i]^2))/sum(memb[,i]^m)
		  	for(j in 1:n) {
	  			dist2cent[j,i] = sqrt((sum((memb[,i]^m)*(x[j,]^2))/sum(memb[,i]^m))-vargeom)
			  }
		  }
	  } else { # For medoid methods, medoids are stored in mobileCenters and fixedCenters 
       med = c(y$mobileCenters,y$fixedCenters)
	     for(i in 1:k) {
	       dist2cent[,i] = x[,med[i]] 
	      }
	  }
	}
	
	  #2. compute membership to centroids
	if (y$method=="KM"||y$method=="KMdd") {
	  minC<-apply(dist2cent,1,which.min)
	  u[,] = 0
	  for(k in 1:length(minC)) u[k,minC[k]] = 1.0
	} else if(y$method=="NC") {
	  d2cm2<-cbind(dist2cent,dnoise)^2
	  for(k in 1:ncol(d2cm2)) {
	    a<-sweep(d2cm2,1,d2cm2[,k],"/")
	    u[,k] = 1/rowSums(a^(-1/(m-1)))
	  }
	  u[d2cm2==0]=1
	} else if(y$method=="HNC"||y$method=="HNCdd") {
	  d2cm<-cbind(dist2cent,dnoise)
	  u[,] = 0
	  minC<-apply(d2cm,1,which.min)
	  for(k in 1:length(minC)) {
	    u[k,minC[k]] = 1.0
	  }
	} else if(y$method=="NCdd") {
	  d2cm<-cbind(dist2cent,dnoise)
	  for(k in 1:ncol(d2cm)) {
	    a<-sweep(d2cm,1,d2cm[,k],"/")
	    u[,k] = 1/rowSums(a^(-1/(m-1)))
	  }
	  u[d2cm==0]=1
	} else if (y$method=="FCM") {
	  d2cm2<-dist2cent^2
	  for(k in 1:ncol(dist2cent)) {
	    a<-sweep(d2cm2,1,d2cm2[,k],"/")
	    u[,k] = 1/rowSums(a^(-1/(m-1)))
	  }
	  u[dist2cent ==0]=1
	} else if (y$method=="FCMdd") {
	  d2cm<-dist2cent
	  for(k in 1:ncol(d2cm)) {
	    a<-sweep(d2cm,1,d2cm[,k],"/")
	    u[,k] = 1/rowSums(a^(-1/(m-1)))
	  }
	  u[dist2cent ==0]=1
	} else if (y$method=="PCM") {
	  for(k in 1:ncol(dist2cent)) u[,k] = 1/(1+((dist2cent[,k]^2)/eta[k])^(1/(m-1)))
	  u[dist2cent==0]=1
	} else if (y$method=="PCMdd") {
	  for(k in 1:ncol(dist2cent)) u[,k] = 1/(1+((dist2cent[,k])/eta[k])^(1/(m-1)))
	  u[dist2cent==0]=1
	} 	 	
	
   #Prepare output
   u = as.data.frame(u)   
   dist2cent = as.data.frame(dist2cent)   
   if(ncol(u)==ncol(y$memb)) names(u) = names(y$memb)
   else {
   	names(u)[1:ncol(y$memb)] = names(y$memb)
   	names(u)[ncol(y$memb)+1] = "N"
   }
   names(dist2cent) = names(y$dist2clusters)   
	rownames(u) = rownames(x)
	rownames(dist2cent) = rownames(x)
	
   res = list(method = y$method, m =y$m, dnoise = y$dnoise, eta= y$eta, memb=u,dist2clusters=dist2cent)
   class(res)<-"vegclass"
	return(res)
		
}