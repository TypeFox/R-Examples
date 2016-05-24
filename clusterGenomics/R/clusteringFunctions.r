
####################################################################################################################
## Author: Gro Nilsen
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the clusterGenomics package
## Reference: "Identifying clusters in genomics data by recursive partitioning", Nilsen et al. (2013, preprint)
####################################################################################################################

# A range of helper functions for application of PART and gap #

### Required by: ##
## part.r
## gap.r


#Compute distance matrix for given distance measure; note: distance between rows!
getDist <- function(X,dist.method,cor.method="pearson"){
	#Calculate distance matrix:
	if(dist.method=="sq.euclidean"){
		dX <- dist(X,method="euclidean")^2
	}else if(dist.method=="cor"){
		#cor computes correlation between columns
		dX <- as.dist(1-cor(t(X),method=cor.method))
	}else{
		dX <- dist(X,method=dist.method)
	}

	return(dX)

}

#Perform hierarchical clustering with given linkage, and apply horizontal cutting of tree into k clusters
doHclust <- function(d,k,linkage){
	cl <- hclust(d,method=linkage)
	lab <- cutree(cl,k)
	return(list(cl=cl,lab=lab))
}

#Perform K-means clustering into k clusters
doKmeans <- function(X,k,nstart){
	cl <- kmeans(X,k,nstart=nstart)                        #distance is always euclidean
  lab <- cl$cluster
	return(list(cl=cl,lab=lab))
}


#Obtain partitions of the data into k=1,...,Kmax clusters
#... gives the clustering method, distance measure, and other parameters to be used in the clustering
findPartition <- function(X,Kmax,dX=NULL,...){

  arg <- as.list(...)
  
	#Calculate distances:
	if(is.null(dX) && arg$cl.method!="kmeans"){
    dX <- getDist(X,dist.method=arg$dist.method,cor.method=arg$cor.method)
	}
	
	#Find the partition for k=1,..,Kmax and return cluster labels stored in a list
  cl.lab <- vector("list",Kmax)
  for(k in 1:Kmax){
    labX <- switch(arg$cl.method,
				  hclust=doHclust(dX,k=k,linkage=arg$linkage)$lab,
				  kmeans=doKmeans(X,k,nstart=arg$nstart)$lab)    #note:kmeans only calculated for euclidean distance!
				  
    cl.lab[[k]] <- labX
  }   
  return(cl.lab)		
}



