
####################################################################################################################
## Author: Gro Nilsen
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the clusterGenomics package
## Reference: "Identifying clusters in genomics data by recursive partitioning", Nilsen et al. (2013, preprint)
####################################################################################################################


### Requires: ##
## clusteringFunctions.r
## gap.r

#Main function for PART:

part <- function(X,Kmax=10,minSize=8,minDist=NULL,cl.lab=NULL,...){

  #Default ... values:
  default.par <- list(q=0.25,Kmax.rec=5,B=100,ref.gen="PC",dist.method="euclidean",cl.method="hclust",linkage="average",cor.method="pearson",nstart=10)
  #Check for user modifications:
  fixed.par <- c(minDist=minDist,minSize=minSize,modifyList(default.par,list(...)))

  
  #Find stopping threshold if minDist is NULL
  if(is.null(minDist)){
    minDist <- get.threshold(X,q=fixed.par$q,fixed.par)
    fixed.par$minDist <- minDist
  }
  
  #Start recursive runs:
  clusters = PartRec(X,Kmax=Kmax,ind=rep(1,nrow(X)),cl.lab=cl.lab,fixed.par)
  
  
  #Check for possible outliers and assign cluster labels
  label <- getPARTlabels(clusters,minSize)
  outliers <- which(label==0)
  if(length(outliers)==0){
    outliers <- NULL
    hatK <- length(unique(label))
  }else{
    hatK <- length(unique(label[-outliers]))
  }
  
  
  return(list(hatK=hatK,lab.hatK=label,outliers=outliers))
}




#The recursive function:

PartRec <- function(X,Kmax,ind,cl.lab=NULL,...){

  fixed.par <- as.list(...)

  #STEP 1: Make sure it is feasible to split X into two clusters each of size >= minSize, otherwise return this cluster
  if(sum(ind)<(2*fixed.par$minSize)){
    return(ind)
  }
  

  #STEP 2: Use a clustering algorithm to partition the objects in X into K=1,..,Kmax clusters: 
	#First Make sure Kmax does not exceed the number of objects in X:
  n <- sum(ind)
	if(n<=Kmax){
    if(fixed.par$cl.method=="hclust"){
      Kmax = n
    }else{
      Kmax = n-1   #kmeans does not work if K=n
    }
  }
  if(is.null(cl.lab)){
	 cl.lab <- findPartition(X=X,Kmax=Kmax,dX=NULL,fixed.par)
	}
  
  
  
  #STEP 3: #Use an objective function (Gap) to decide on the optimal number of clusters, hatK, for this set X
  gap.res <- gap(X=X,Kmax=Kmax,cl.lab=cl.lab,B=fixed.par$B,ref.gen=fixed.par$ref.gen,fixed.par=fixed.par)
  hatK <- gap.res$hatK
  lab.hatK <- gap.res$lab.hatK
  
  
  #STEP 4.
  
  #In case hatK > 1, make sure at least two of them are >= minSize:
  if(sum(table(lab.hatK) >= fixed.par$minSize)<2){
    hatK <- 1  
  } 
  
  
  #a) hatK == 1:
  if(hatK==1 && length(cl.lab)==1){
    #Special case if minSize=1 and only 2 obs in X; kmeans cannot return 2 clusters (see findPartition) and Kmax is therefore set to 1 -> cannot divide into two tentative clusters
    return(ind)
  }
  
  if(hatK==1){
    #Divide set into two tentative clusters:
    obs1 <- cl.lab[[2]]==1
    obs2 <- cl.lab[[2]]==2
    T1 <- X[obs1,,drop=FALSE]   #drop is necessary in case obs1 is of length 1!
    T2 <- X[obs2,,drop=FALSE]
    
    
    # Check if stopping threshold has been reached:
    #Use hierarchical clustering to determine the distance between T1 and T2:
    hc.res <- doHclust(getDist(X,dist.method=fixed.par$dist.method,cor.method=fixed.par$cor.method),k=1,linkage=fixed.par$linkage)$cl    #(k is irrelevant here, only specified because doHclust needs it)
	  T.height <- max(hc.res$height) 
     
  
    if(T.height > fixed.par$minDist){
      #Create new index-vectors corresponding to the two tentative clusters:
      t1.ind <- ind
      t2.ind <- ind
      t1.ind[which(ind==1)[!obs1]] <- 0
      t2.ind[which(ind==1)[!obs2]] <- 0
      #Tentative runs:
      t1 = PartRec(X=T1,Kmax=fixed.par$Kmax.rec,ind=t1.ind,cl.lab=NULL,fixed.par)
      t2 = PartRec(X=T2,Kmax=fixed.par$Kmax.rec,ind=t2.ind,cl.lab=NULL,fixed.par)
      #If no clusters are found in recursive runs the return value will be a vector, otherwise it will be a matrix:
      if(!is.matrix(t1) && !is.matrix(t2)){
        return(ind)
      }else{
        return(cbind(t1,t2))
      }
    }else{
      #Return the current set if stopping criterion has been reached
      return(ind)
    }
  
  }else{
    #b) hatK > 1
    
    res <- matrix(NA,nrow=length(ind),ncol=0)
    for(k in 1:hatK){
      #Pick out subset in this cluster
      obs = cl.lab[[hatK]]==k
      S = X[obs,,drop=FALSE]
      ind.S = ind
      ind.S[which(ind==1)[!obs]] <- 0 
      res = cbind(res,PartRec(X=S,Kmax=fixed.par$Kmax.rec,ind=ind.S,cl.lab=NULL,fixed.par))
    }
    return(res)

  }

}



#Help-functions only used by part:

#Find a stopping threshold given the percentage of heights to be used in the dendrogram
get.threshold <- function(X,q,...){
	arg <- as.list(...)
  
	#Get distance matrix
	dX <- getDist(X,dist.method=arg$dist.method,cor.method=arg$cor.method)
	
	#Get hierarchical clustering result:
	cl <- hclust(dX,method=arg$linkage)

	#The total set of cluster heights in dendrogram:
	h <- cl$height
	
  use.h <- quantile(h,probs=1-q)
  
  return(use.h)

}



getPARTlabels <- function(clusters,minSize){
  if(!is.matrix(clusters)){
    clusters <- as.matrix(clusters)
  }
  label <- rep(NA,nrow(clusters))
  ncl <- ncol(clusters)
  id=1
  for(j in 1:ncl){
    if(sum(clusters[,j])<minSize){
      label[clusters[,j]==1] <- 0
    }else{
      label[clusters[,j]==1] <- id
      id=id+1
    }                            
  
  }
  return(label)
}





