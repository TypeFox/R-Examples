defuzzify<-function(object, method="max", alpha=0.5,na.rm=FALSE) {
    METHODS <- c("max", "cut")
    method <- match.arg(method, METHODS)	
	if(inherits(object,"vegclust")|| inherits(object,"vegclass")) memb = object$memb
	else memb = object
	
	memb = as.data.frame(memb)
	clnames = names(memb)
	cluster = vector("character",nrow(memb))
	u = matrix(0,nrow(memb), ncol(memb))
	cluster = rep(NA,nrow(memb))
	
	a<-function(v) {
		return(paste(clnames[which(v==1)],collapse="+"))
	}
  if(method=="max") {
  	  vmax<-apply(memb,1,max)
   	u=ifelse(memb==vmax,1,0)
  } else if (method=="cut") {
   	u[memb>alpha]= 1
  }
  cluster = as.character(apply(u,1,a))
  cluster[cluster==""]= NA

	u = as.data.frame(u)
	row.names(u) = row.names(memb)
	names(u) = names(memb)
	if(na.rm) {
		sel = !is.na(cluster)
		cluster =cluster[sel]
		u = subset(u,sel)
	}
  names(cluster) = row.names(u)
	return(list(memb=u,cluster=cluster))
}