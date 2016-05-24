index.H<-function(x,clall,d=NULL,centrotypes="centroids")
{
wgss<-function(x,cl,d,centroids)
{
        n <- length(cl)
        k <- max(cl)
  if(is.null(dim(x))){
    dim(x)<-c(length(x),1)
  }
	centers<-matrix(nrow=k,ncol=ncol(x))
  for(i in 1:k){
    if(centrotypes=="centroids"){
      if(ncol(x)==1){
      centers[i,]<-mean(x[cl==i,])
      }
      else{
        if(is.null(dim(x[cl==i,]))){
          centers[i,]<-x[cl==i]
        }
        else{
          centers[i,]<-apply(x[cl==i,],2,mean)
        }
      }
    }
    else{
      centers[i,]<-.medoid(x[cl==i,],d[cl==i,cl==i])
   }
  }
     withins <- rep(0, k)
	   x <- (x - centers[cl,])^2
        for(i in 1:k){
           withins[i] <- sum(x[cl==i,])
	   }
        sum(withins)
}
   	 if(sum(c("centroids","medoids")==centrotypes)==0)
      stop("Wrong centrotypes argument")
   	 if("medoids"==centrotypes && is.null(d))
      stop("For argument centrotypes = 'medoids' d cannot be null")
     if(!is.null(d)){
      if(!is.matrix(d)){
        d<-as.matrix(d)
      }
     row.names(d)<-row.names(x)
     }
        if(is.null(dim(x))){
          dim(x)<-c(length(x),1)
        }
        n <- nrow(x)
        g <- max(clall[,1])
        (wgss(x,clall[,1],d,centrotypes)/wgss(x,clall[,2],d,centrotypes)-1)*(n-g+1)
}

