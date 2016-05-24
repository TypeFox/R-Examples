#' Information Gain 
#' 
#' Compute the cost function of a node of samples in a tree 
#'  
#' @param y Output Features for the samples of the node 
#' @param V_inv Covariance matrix of Output Feature matrix
#' @param Command 1 for RF and 2 for MRF depending on the method
#' @return cost or Entropy of samples in a node of a tree 
#' @details
#' In multivariate trees(MRF) node cost is measured as the sum of squares of the Mahalanobis 
#' distance to capture the correlations in the data where in univariate trees node cost is measured as the Euclidean distance. 
#' Mahalanobis Distance captures the distance of the sample point from the mean of the node along the principal component axes.
#' @examples
#' y=matrix(runif(10*2),10,2)
#' V_inv=stats::cov(y)
#' Command=2
#' #Command=2 for MRF and 1 for RF
#' #This function calculates information gain of a node
#' Cost=Multi_D_mod(y,V_inv,Command)
#' @export
#' @references
#' De Maesschalck, Roy, Delphine Jouan-Rimbaud, and Desire L. Massart. "The mahalanobis distance." Chemometrics and intelligent laboratory systems 50.1 (2000): 1-18.
Multi_D_mod  <- function(y,V_inv,Command){
  NN=nrow(y)
  ybar5=colSums (y, na.rm = FALSE, dims = 1)/nrow(y)
  ybar2=matrix(ybar5,nrow=1,ncol=ncol(y))
  
  if (Command==2){ #using VMRF
    yhat=y-kronecker(matrix(1,nrow(y),1),ybar2)
    
    D = sum(diag(yhat %*% V_inv %*% t(yhat)))
    
  }else if(Command==1){  #using rf
    ybar=sum(y)/nrow(y)
    D=sum((y-ybar)^2)
  }
  return(D)
}