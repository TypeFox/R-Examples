.medoid<-function(x,d)
{
  minj<-0
  minsumdist<-sum(d)
  if(is.null(dim(x)) && is.null(dim(d))){
     dim(x)<-c(1,length(x))
     x
  }
  else{
    if(is.null(dim(d))){
          dim(d)<-c(1,1)
    }
    if(is.null(dim(x))){
      dim(x)<-c(length(x),1)
    }
    for(j in 1:nrow(d)){
      if (sum(d[j,])<=minsumdist){
         #minj<-row.names(d)[j]
         minj<-j
         minsumdist<-sum(d[j,])
      }
    }
    resul<-as.matrix(x[minj,])  
    resul
  }
}


index.G1<-function(x,cl,d=NULL,centrotypes="centroids")
{
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
     
        n <- length(cl)
        k <- max(cl)
        if(is.null(dim(x))){
          dim(x)<-c(length(x),1)
        }
	   centers<-matrix(nrow=k,ncol=ncol(x))
      for(i in 1:k)
	   {
    x.k = x[cl==i,]
    if(centrotypes=="centroids"){
      if(ncol(x)==1){
        centers[i,]<-mean(x.k)
      }
      else{
        if (is.vector(x.k)){
          centers[i,]<-x.k
        }
        else{
          centers[i,]<-apply(x.k,2,mean)
        }
      }
    }
    else{
        centers[i,]<-.medoid(x[cl==i,],d[cl==i,cl==i])
        #print(apply(x[cl==i,],2,mean))
        #print(centers[i,])
       }
	   }
        if (centrotypes=="centroids"){
          allmean <- apply(x,2,mean)
        }
        else{
          # print(apply(x,2,mean))
          allmean<-.medoid(x,d)
          #print(allmean)
        }
        dmean <- sweep(x,2,allmean,"-")
        allmeandist <- sum(dmean^2)
        withins <- rep(0, k)
	   x <- (x - centers[cl,])^2
        for(i in 1:k){
           withins[i] <- sum(x[cl==i,])
	   }
        wgss <- sum(withins)
        bgss <- allmeandist - wgss
        (bgss/(k-1))/(wgss/(n-k))

}

