#' CrossClustering: a Partial Clustering Algorithm with Automatic Estimation of the Number of Clusters and Identification of Outliers
#'
#' This function performs the CrossClustering algorithm. This method combines the Ward's minimum variance and Complete Linkage algorithms, providing automatic estimation of a suitable number of clusters and identification of outlier elements.
#'
#' @param d a dissimilarity structure as produced by the function \code{dist}
#' @param k.w.min minimum number of clusters for the Ward's minimum variance method. By default is set equal 2
#' @param k.w.max maximum number of clusters for the Ward's minimum variance method (see details)
#' @param k.c.max maximum number of clusters for the Complete-linkage method. It can not be equal or greater than the number of elements to cluster (see details)
#' @param out logical. If \code{TRUE} (default) outliers must be searched (see details)
#' @return A list of objects describing characteristics of the partitioning as follows:
#' \item{Optimal.cluster}{number of clusters}
#' \item{Cluster.list}{a list of clusters; each element of this lists contains the indices of the elemenents belonging to the cluster}
#' \item{Silhouette}{the average silhouette witdh over all the clusters}
#' \item{n.total}{total number of input elements}
#' \item{n.clustered}{number of input elements that have actually been clustered}
#'  
#'  @details See cited document for more details.
#'  @examples
#'  
#' ### Generate simulated data
#' toy <- matrix(NA, nrow=10, ncol=7)
#' colnames(toy) <- paste("Sample", 1:ncol(toy), sep="")
#' rownames(toy) <- paste("Gene", 1:nrow(toy), sep="")
#' set.seed(123)
#' toy[,1:2] <- rnorm(n=nrow(toy)*2, mean=10, sd=0.1)
#' toy[,3:4] <- rnorm(n=nrow(toy)*2, mean=20, sd=0.1)
#' toy[,5:6] <- rnorm(n=nrow(toy)*2, mean=5, sd=0.1)
#' toy[,7] <- runif(n=nrow(toy), min=0, max=1)
#' 
#' ### toy is transposed as we want to cluster samples (columns of the original matrix)
#' d <- dist(t(toy), method="euclidean")
#' 
#' ### Run CrossClustering
#' toyres <- CrossClustering(d, k.w.min=2, k.w.max=5, k.c.max=6, out=TRUE)
#'  
#'  @author
#'  Paola Tellaroli, \email{paola.tellaroli@unipd.it}; Michele Donato, \email{michele.donato@wayne.edu}
#'  
#'  @references
#'  Tellaroli, P., Bazzi, M., Donato, M., Brazzale, A. R., Draghici, S. (2016), Cross Clustering: a Partial Clustering Algorithm with Automatic Estimation of the Number of Clusters, Plos ONE (In Press)

CrossClustering<-function(d,k.w.min=2,k.w.max,k.c.max,out=TRUE)
{
  #require(cluster)
  n <- (1+sqrt(1+8*length(d)))/2
  beta.clu.ward<- hclust(d, method="ward.D")
  beta.clu.complete<-hclust(d, method="complete")
  grid<-as.matrix(expand.grid(k.w.min:k.w.max,k.w.min:k.c.max))
  if (out==T) 
    grid<-grid[grid[,2]>grid[,1],] 
  else grid<-grid[grid[,2]>=grid[,1],]
  grid<-cbind(grid,0)
  colnames(grid)<-c("Ward","Complete","N. classified")
  n.clu <- NULL
  for(i in 1:dim(grid)[1]){
    n.clu[i] <- proportion.function(grid[i,],  beta.clu.ward=beta.clu.ward,
                                        beta.clu.complete=beta.clu.complete)
  }
  grid[,3] <- n.clu
  grid.star<-which(grid== max(grid[,3], na.rm = FALSE), arr.ind = TRUE)[,1] 
  k.star<-rbind(grid[grid.star,1:2]) 
  if(is.null(dim(k.star))){
    cluster.list <- proportion.function(k.star, beta.clu.ward=beta.clu.ward, 
                                            beta.clu.complete=beta.clu.complete, return.list=T)
    clustz <- sapply(1:n, geneinlista, cluster.list$beta.list)
  }else{
    cluster.list <- apply(k.star, 1, proportion.function, beta.clu.ward=beta.clu.ward, beta.clu.complete=beta.clu.complete, return.list=T)
    clustz <- sapply(cluster.list, function(lasim) sapply(1:n, geneinlista, lista=lasim$beta.list)) 
    
  }
  clustz[clustz=="integer(0)"]=0
  if(is.null(dim(clustz))){
    clustz <- matrix(clustz, ncol=1)
  }
  Sil <- list()
  for (c in 1:ncol(clustz)){
    Sil[c] <- mean(silhouette(as.numeric(clustz[,c]), dist=d)[,3])
  }
  
  if(is.null(dim(k.star))){
    k.star.star <- k.star[which.max(Sil)]
  }else{
    k.star.star <- k.star[which.max(Sil),]
  }
  Cluster.list <- cluster.list[[which.max(Sil)]]$beta.list
  n.clustered <-length(unlist(Cluster.list))
  return(list("Optimal.cluster"=length(cluster.list[[which.max(Sil)]]$beta.list), 
              "Cluster.list"=Cluster.list, 
 #             "A.star"=cluster.list[[which.max(Sil)]]$A.star,
              "Silhouette"=max(unlist(Sil), na.rm = FALSE),
              "n.total"=n,
              "n.clustered"=n.clustered)
  )
}
