#' Make test data
#' 
#' \code{makeData} Function for making test data
#' 
#' @param p Number of features (or variables)
#' @param n Number of observations
#' @param k Number of clusters
#' @param sigma Variance
#' @param missing Desired missingness percentage
#' @param seed (Optional) Seed (default seed is 12345)
#' 
#' @export
#' 
#' @examples
#' p <- 2
#' n <- 100
#' k <- 3
#' sigma <- 0.25
#' missing <- 0.05
#' 
#' X <- makeData(p,n,k,sigma,missing)$Orig
#' 
#' @author Jocelyn T. Chi
makeData <- function(p,n,k,sigma,missing,seed=12345){
  if(p <= 0 | n <= 0 | k <= 0){
    return(cat('Please select positive values for p, n, and k.'))
  }
  if (missing < 0 | missing > 1){
    return(cat('Please select a missingness percentage between 0 and 1.'))
  }
  
  # Make complete data
  set.seed(seed)
  M <- matrix(rnorm(k*p),k,p)
  assignment <- sample(1:k,n,replace=TRUE)
  X <- M[assignment,] + sigma*matrix(rnorm(n*p),n,p)
  
  # Make missing data
  X_missing <- X
  missing_ix <- sample(1:(n*p),(n*p*missing),replace=TRUE)
  X_missing[missing_ix] <- NA
  
  return(list(Orig=X,Missing=X_missing,truth=assignment))  
}

#' Function for assigning clusters to rows in a matrix
#' 
#' \code{assign_clustpp} Function for assigning clusters to rows in a matrix
#' 
#' @param X Data matrix containing missing entries whose rows are observations and columns are features
#' @param init_centers Centers for initializing k-means
#' @param kmpp_flag (Optional) Indicator for whether or not to initialize with k-means++
#' @param max_iter (Optional) Maximum number of iterations
#' 
#' @export
#' 
#' @examples
#' p <- 2
#' n <- 100
#' k <- 3
#' sigma <- 0.25
#' missing <- 0.05
#' Data <- makeData(p,n,k,sigma,missing)
#' X <- Data$Missing
#' Orig <- Data$Orig
#' 
#' clusts <- assign_clustpp(Orig, k)
#' 
#' @author Jocelyn T. Chi
assign_clustpp <- function(X,init_centers,kmpp_flag=TRUE,max_iter=20){
  res <- kmeans(X, init_centers)
  clusts <- res$cluster
  obj <- res$totss
  fit <- 1-(sum(res$withinss)/res$totss)
  centers <- res$centers
  if (kmpp_flag == TRUE) {
    ## Try to find a better assignment
    for (iter in 1:max_iter) {
      centers_kmpp <- kmpp(X,length(res$size))
      sol <- kmeans(X, centers_kmpp)
      if (sol$totss < obj) {
        obj <- sol$totss
        clusts <- sol$cluster
        fit <- 1-(sum(sol$withinss)/sol$toss)
        centers <- sol$centers
        break
      }
    }
  }
  return(list(clusts=clusts,obj=obj,centers=centers,fit=fit))
}

#' Function for finding indices of missing data in a matrix
#' 
#' \code{findMissing} Function for finding indices of missing data in a matrix
#' 
#' @param X Data matrix containing missing entries whose rows are observations and columns are features
#' 
#' @return A numeric vector containing indices of the missing entries in X
#' 
#' @export
#' 
#' @examples
#' p <- 2
#' n <- 100
#' k <- 3
#' sigma <- 0.25
#' missing <- 0.05
#' Data <- makeData(p,n,k,sigma,missing)
#' X <- Data$Missing
#' missing <- findMissing(X)
#' 
#' @author Jocelyn T. Chi
findMissing <- function(X){
  missing_all <- which(is.na(X))
  return(missing_all)
}

#' Function for initial imputation for k-means
#' 
#' \code{initialImpute} Initial imputation for k-means
#' 
#' @param X Data matrix containing missing entries whose rows are observations and columns are features
#' 
#' @return A data matrix containing no missing entries
#' 
#' @export
#' 
#' @examples
#' p <- 2
#' n <- 100
#' k <- 3
#' sigma <- 0.25
#' missing <- 0.05
#' Data <- makeData(p,n,k,sigma,missing)
#' X <- Data$Missing
#' X_copy <- initialImpute(X)
#' 
#' @author Jocelyn T. Chi
initialImpute <- function(X){
  avg <- mean(X,na.rm=TRUE)
  X[which(is.na(X))] <- avg
  return(X)
}

#' Function for performing k-POD
#' 
#' \code{kpod} Function for performing k-POD, a method for k-means clustering on partially observed data
#' 
#' @param X Data matrix containing missing entries whose rows are observations and columns are features
#' @param k Number of clusters
#' @param kmpp_flag (Optional) Indicator for whether or not to initialize with k-means++
#' @param maxiter (Optional) Maximum number of iterations
#' 
#' @return cluster: Clustering assignment obtained with k-POD
#' @return cluster_list: List containing clustering assignments obtained in each iteration
#' @return obj_vals: List containing the k-means objective function in each iteration
#' @return fit: Fit of clustering assignment obtained with k-POD (calculated as 1-(total withinss/totss))
#' @return fit_list: List containing fit of clustering assignment obtained in each iteration
#' 
#' @export
#' 
#' @import clues
#' 
#' @examples
#' p <- 5
#' n <- 200
#' k <- 3
#' sigma <- 0.15
#' missing <- 0.20
#' Data <- makeData(p,n,k,sigma,missing)
#' X <- Data$Missing
#' Orig <- Data$Orig
#' truth <- Data$truth
#' 
#' kpod_result <- kpod(X,k)
#' kpodclusters <- kpod_result$cluster
#' 
#' @author Jocelyn T. Chi
#' 
kpod <- function(X,k,kmpp_flag=TRUE,maxiter=100){
  
  n <- nrow(X)
  p <- ncol(X)
  
  cluster_vals <- vector(mode="list",length=maxiter)
  obj_vals <- double(maxiter)
  fit <- double(maxiter)
  
  missing <- findMissing(X)
  
  # Run first iteration
  X_copy <- initialImpute(X)
  
  ## Use kmpp to select initial centers
  init_centers <- kmpp(X_copy, k)
  temp <- kmeans(X_copy,init_centers)
  clusts <- temp$cluster
  centers <- temp$centers
  fit[1] <- 1-(sum(temp$withinss)/temp$totss)
  
  # Make cluster matrix
  clustMat <- centers[clusts,]
  
  # Update with centers
  X_copy[missing] <- clustMat[missing]
  
  #obj_vals[1] <- temp$obj
  obj_vals[1] <- sum((X[-missing]-clustMat[-missing])^2)
  cluster_vals[[1]] <- clusts
  
  # Run remaining iterations
  for (i in 2:maxiter){
    temp <- assign_clustpp(X_copy,centers,kmpp_flag)
    clusts <- temp$clusts
    centers <- temp$centers
    fit[i] <- temp$fit
    
    # Impute clusters
    clustMat <- centers[clusts,]
    X_copy[missing] <- clustMat[missing]
    
    obj_vals[i] <- sum((X[-missing]-clustMat[-missing])**2)
    cluster_vals[[i]] <- clusts
    
    if (all(cluster_vals[[i]] == cluster_vals[[i-1]])){
      noquote('Clusters have converged.')
      return(list(cluster=clusts,cluster_list=cluster_vals[1:i],obj_vals=obj_vals[1:i],fit=fit[i],fit_list=fit[1:i]))
      break
    }
  }
  return(list(cluster=clusts,cluster_list=cluster_vals[1:i],obj_vals=obj_vals[1:i],fit=fit[i],fit_list=fit[1:i]))
}

#' k-means++
#'
#' \code{kmpp} Computes initial centroids via kmeans++
#'
#' @param X Data matrix whose rows are observations and columns are features
#' @param k Number of clusters.
#' 
#' @return A data matrix whose rows contain initial centroids for the k clusters
#' 
#' @export
#' 
#' @examples
#' n <- 10
#' p <- 2
#' X <- matrix(rnorm(n*p),n,p)
#' k <- 3
#' kmpp(X,k)
#' 
kmpp <- function(X, k) {
  n <- nrow(X)
  p <- ncol(X)
  C <- integer(k)
  C[1] <- sample(1:n, 1)
  for (i in 2:k) {
    S <- matrix(NA,n,i-1)
    for (j in 1:(i-1)) {
      S[,j] <- apply(X -
                       matrix(X[C[j],],n,p,byrow=TRUE),1,FUN=function(x)
                       {norm(as.matrix(x),'f')**2})
    }
    D <- apply(S,1,min)
    pr <- D/sum(D)
    C[i] <- sample(1:n, 1, prob = pr)
  }
  return(X[C,])
}