#' Nearest neighbor
#' 
#' An auxiliary  function to find the nearest neighbors from the distance matrix
#' 
#' 
#' @param model nearest neighbor ranking model 
#' @param k number of the nearest neighbors to consider
#' @return a vector of length of \code{model}
#' @details This function is applied to find the nearest neighbors in the distance matrix.
kNN <- function(model, k) {
  w <- rep(0, length(model))
  w[order(model)[1:k]] <- 1
  colnames(w) <- colnames(model)
  w
}


#' Predicting rankings using the nearest neighbor algorithm
#' 
#' This function makes prediction of rankings based on the nearest neighbor
#' @param train.x is matrix of numeric attributes in training sample
#' @param test.x is a vector of new numeric attributes for which to predict rankings
#' @param y is matrix of training rankings
#' @param k is the number of the nearest neighbors to consider (default k=3)
#' @param n is a parameter of 'memory' of how fast the past rankings gets forgotten. (see details of \link{time_weights}). By default, \code{n=1} which means that a label ranking problem does not have timing effect. 
#' @return a vector of predicted ranking for attribute \code{test.x}
#' @details A function predicts the rankings based on the euclidean distance between train and test attributes. 
#' 
#' @importFrom pdist pdist
#' 
#' @examples 
#' train.x <- lr.num[1:16,]
#' test.x <- lr.num[17,]
#' ranking <- nn_rank(train.x, y, test.x, n=1,k=3)
#' @export
nn_rank <- function(train.x, y, test.x, n=1, k=3) {
  w <- time_weights(nrow(as.matrix(train.x)),n)
  model <- w * as.matrix(pdist(test.x, train.x))
  nn <- kNN(model,k)
  rank(apply(y[which(nn == 1), drop = F, ], 2, mean))
}