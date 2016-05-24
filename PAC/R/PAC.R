#' @importFrom stats kmeans mahalanobis 
#' @importFrom utils combn
#' @import Rcpp
#' @useDynLib PAC
NULL

#' Finds N Leaf centers in the data
#' 
#' @param  data     a n x p data matrix
#' @param  N        number of leaves centers
#' @param  method   partition method, either "dsp(discrepancy based partition)",  or "bsp(bayesian sequantial partition)"
#' @return leafctr  N leaves centers
#' @export
#' @import Rcpp
BSPLeaveCenter <- function(data, N = 40, method = "dsp") {
  leafctr = BSPLeaveCenterCpp(data, N, method) # cpp function
  return(leafctr)
} 


#' PAC (Partition Assisted Clustering)
#' 
#' @param  data     a n x p data matrix
#' @param  K        number of final clusters in the output
#' @param  maxlevel the maximum level of the partition
#' @param  method   partition method, either "dsp(discrepancy based partition)",  or "bsp(bayesian sequantial partition)"
#' @param  max.iter maximum iteration for the kmeans step
#' @return y        cluter labels for the input
#' @export
PAC <- function(data, K, maxlevel = 40, method = "dsp", max.iter = 50) {
  if (typeof(data) != "matrix")
    data = as.matrix(data)
  
  dim = ncol(data)
  clusterDist <- function(C) {
    x_idx = IDX == C[1]
    y_idx = IDX == C[2]
    nx = sum(x_idx);
    ny = sum(y_idx);
   
    if (dim < 2) {
      xbar = mean(data[x_idx,])
    } else if(nx > 1) {
      xbar = colMeans(data[x_idx,]);
    } else {
      xbar = data[x_idx,];
    }
    
    if (dim < 2) {
      ybar = mean(data[y_idx,])
    } else if(ny > 1) {
      ybar = colMeans(data[y_idx,]);
    } else {
      ybar = data[y_idx,];
    }
    
    if (dim < 2) {
      d = sqrt((xbar - ybar)^2)
    } else if (nx <= ncol(data) || ny <= ncol(data)) {
      d = sqrt(sum((xbar - ybar)^2))
    } else {
      Sx = var(data[x_idx,])
      Sy = var(data[y_idx,])
      if (kappa(Sx) > 1e5 || kappa(Sy) > 1e5) {
        d = sqrt(sum((xbar - ybar)^2))
      } else {
        d = min(mahalanobis(xbar, ybar, Sx), mahalanobis(xbar, ybar, Sy))
      }
    }
    return(d)
  }
  
  leafcenters = BSPLeaveCenter(data, maxlevel, method)
  
  maxC = nrow(leafcenters)
  print(maxC)
  print("Initial Clustering...")
  options(warn=-1)
  fit = kmeans(data, leafcenters, algorithm = "Lloyd", iter.max = max.iter, trace = TRUE)
  options(warn=0)
  cluster.size = fit$size
  print("Merging...")
  IDX = fit$cluster
  label = unique(IDX)
  numcluster = length(label)
  clusts = t(combn(label,2))
  start = proc.time()
  D = cbind(clusts, apply(clusts, 1, clusterDist))
  
  while (numcluster > K) {
    min_idx = which.min(D[,3]);
    min_i = D[min_idx,1];
    min_j = D[min_idx,2];
    # remove cluster min_j
    IDX[IDX == min_j] = min_i;
    # remove the cluster indexed by min_j
    rm_idx = D[, 2] == min_j | D[, 1] == min_j
    D[rm_idx,] <- matrix();
    D = D[complete.cases(D),]
    rc_idx = which(D[,1] == min_i | D[,2] == min_i)
    for (k in 1:length(rc_idx)) {
      D[rc_idx[k], 3] = clusterDist(c(D[rc_idx[k], 1], D[rc_idx[k], 2]))
    }
    numcluster = length(unique(IDX))
  }
  y = IDX;
  newlevel = unique(IDX)
  for (i in 1:length(newlevel)) {
    idx = IDX == newlevel[i];
    y[idx] = i;
  }
  
  return(y)
}



