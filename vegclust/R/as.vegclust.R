as.vegclust <-
function(x,y, method="KM", m=1.0, dnoise=NULL, eta=NULL) {
  METHODS <- c("KM", "FCM", "PCM","NC","HNC" ,"KMdd","NCdd", "HNCdd", "FCMdd", "PCMdd")
  method <- match.arg(method, METHODS)
  
   dist2onecluster<-function(x,object) {
		#x is an (euclidean) distance matrix
		x = as.matrix(x)
		n = nrow(x)
		if (length(object)!=n) stop("Length of object vector must be equal to the number of sites in x")
		vargeom = (rep(1,n) %*% (x^2) %*% rep(1,n))/(2*(n^2))
		return(sqrt((sum(object^2)/n)-vargeom))
	}
   dist2clusters<-function(x,cluster,object) {
	   n = nrow(as.matrix(x))
		if (length(cluster)!=n) 
            stop("Length of cluster vector must be equal to the number of sites in x")
		cluster = as.factor(cluster)
		k = length(levels(cluster))
		d = vector("numeric",k)
		for(i in 1:k) {
			sel = (cluster==levels(cluster)[i])
		  sel[is.na(sel)]=FALSE
			d[i] = dist2onecluster(as.dist(as.matrix(x)[sel,sel]),object[sel])
		}	
		names(d) = levels(cluster)
		return(d)
   }	
	
   if(inherits(x,"dist")) {
   	mode="distance"
   	x = as.matrix(x)
   	sitenames = rownames(x)
   } else {
   	mode="raw"
   	sitenames = rownames(x)
   	varnames = names(x)
   	x = as.matrix(x)
   }
   if(is.vector(y)) {
   	cluster = y
   	cln =levels(as.factor(cluster))
    u = as.memb(cluster)
   	rownames(u) = sitenames
    colnames(u) = cln
   	k = length(cln)
   	if(method=="NC"||method=="NCdd"|| method=="HNC"||method=="HNCdd") {
       u = cbind(u, rep(0,nrow(u)))
       colnames(u)[k+1] = "N"
   	}
    u = as.data.frame(u)
   } else if(is.matrix(y) || is.data.frame(y)) {
   	u = as.data.frame(y)
   	cln = names(u)
   	k = length(cln)
   	if(method=="NC"||method=="NCdd"|| method=="HNC"||method=="HNCdd") {
   	  k = k-1
   	  cln = cln[1:k]
   	}
   }
   n = nrow(x)
   
   dist2cent = matrix(0,nrow=n,ncol=k) 
      
   if(mode=="distance") {
   	for(j in 1:n) {
   		dist2cent[j,] = dist2clusters(x,cluster,x[j,])
   	}
   	centers=NULL
   } else if (mode=="raw"){
   	cm = clustcentroid(x,u[,1:k])
   	colnames(cm) = varnames
   	rownames(cm) = cln
   	centers = as.data.frame(cm)
   	for(i in 1:k) {
   		dist2cent[,i] = sqrt(rowSums(sweep(x,2,cm[i,],"-")^2))
   	}   
   }   
   dist2cent = as.data.frame(dist2cent)  
   names(dist2cent) = cln
   rownames(dist2cent) = sitenames
	
   size = colSums(u[,1:k])
   withinss = colSums((dist2cent^2)*u[,1:k])
   functional = sum(withinss)
	
   res = list(mode = mode, method=method, m = m, dnoise = dnoise, eta = eta, memb=u,mobileCenters=centers, fixedCenters=NULL, dist2clusters=dist2cent, withinss = withinss, size=size, functional=functional)
   class(res)<-"vegclust"
	return(res)
}

