#' Information Gain 
#' 
#' Compute the cost function of a tree node 
#'  
#' @param y Output Features for the samples of the node 
#' @param Cov_Y Covariance matrix of Output Response matrix for MRF(Input 0 for RF)
#' @param Command 1 for univariate Regression Tree (corresponding to RF) and 2 for Multivariate Regression Tree (corresponding to MRF)
#' @return cost or entropy of samples in a node of a tree 
#' @details
#' In multivariate trees (MRF) node cost is measured as the sum of squares of the Mahalanobis distance to capture the correlations in 
#' the data whereas in univariate trees node cost is measured as the sum of Euclidean distance square. Mahalanobis Distance captures 
#' the distance of the sample point from the mean of the node along the principal component axes.
#' @examples
#' y=matrix(runif(10*2),10,2)
#' Cov_Y=stats::cov(y)
#' Command=2
#' #Command=2 for MRF and 1 for RF
#' #This function calculates information gain of a node
#' Cost=Multi_D_mod(y,Cov_Y,Command)
#' @export
#' @references
#' Segal, Mark, and Yuanyuan Xiao. "Multivariate random forests." 
#' Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery 1.1 (2011): 80-87.
Multi_D_mod  <- function(y,Cov_Y,Command){
  NN=nrow(y)
  ybar5=colSums (y, na.rm = FALSE, dims = 1)/nrow(y)
  ybar2=matrix(ybar5,nrow=1,ncol=ncol(y))
  
  if (Command==2){ #using VMRF
    yhat=y-kronecker(matrix(1,nrow(y),1),ybar2)
    
    D = sum(diag(yhat %*% solve(Cov_Y) %*% t(yhat)))
    
  }else if(Command==1){  #using rf
    ybar=sum(y)/nrow(y)
    D=sum((y-ybar)^2)
  }
  return(D)
}