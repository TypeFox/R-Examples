clustvar <-
function(x,cluster=NULL, defuzzify=FALSE,...) {
  "clustvar1"<-function(x) {
	if(inherits(x,"dist")) {
		x = as.matrix(x)
		n = nrow(x)
		return((rep(1,n) %*% (x^2) %*% rep(1,n))/(2*(n^2)))
	}
	else {
		x = as.matrix(x)
		n = nrow(x)
		return(sum(diag(var(x)))*((n-1)/n))
	}
  }
   if(inherits(x,"vegclust")||inherits(x,"vegclass")){
     if(defuzzify) memb = defuzzify(x, ...)$memb
     else memb = x$memb
     k = length(names(x$dist2clusters))
		 v = vector("numeric",k)
   	 for(i in 1:k) {
   	   if(x$method=="KMdd"||x$method=="FCMdd"||x$method=="NCdd"||x$method=="PCMdd"||x$method=="HNCdd") {
   	     v[i] = sum((x$dist2clusters[,i])*(memb[,i]^x$m)) #Distances to medoid weighted by membership
   	   } else {
   	     v[i] = sum((x$dist2clusters[,i]^2)*(memb[,i]^x$m)) #Squared distances to centroid weighted by membership
   	   }
   		 v[i] = v[i]/sum(memb[,i]^x$m)
   	 }
		names(v) = names(x$dist2clusters)
		return(v)		   	
   } else if(is.null(cluster)) {
		return(clustvar1(x))		
	} else{
	   n = nrow(as.matrix(x))
		if (length(cluster)!=n) 
            stop("Length of cluster vector must be equal to the number of sites in x")
		cluster = as.factor(cluster)
		k = length(levels(cluster))
		v = vector("numeric",k)
		for(i in 1:k) {
			sel = (cluster==levels(cluster)[i])
			if(inherits(x,"dist")) {
			   v[i] = clustvar1(as.dist(as.matrix(x)[sel,sel]))
			} else {
			   v[i] = clustvar1(x[sel,])
			}
		}	
		names(v) = levels(cluster)
		return(v)		
	}
}

