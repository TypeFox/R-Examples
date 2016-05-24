#' Model of Single tree of Random Forest or Multivariate Random Forest
#' 
#' Build the model of a tree of RF or MRF using the training samples, which is used for the prediction of testing samples  
#'  
#' @param X Input Training matrix of M x N, M is the number of training samples and N is the number of features
#' @param Y Output Training response of M x T, M is the number of samples and T is number of ouput Features(Response)
#' @param mtree Number of randomly selected features used for each split
#' @param min_leaf Minimum number of samples in the leaf node 
#' @param V_inv Covariance matrix of Output Feature matrix
#' @param Command 1 for RF and 2 for MRF depending on the method
#' @return Model of single tree of the forest(RF or MRF) 
#' @details
#' The Model contains the list of lists, where each list contains either the splitting criteria of a split in
#' a node or, output feature matrix of the samples of a node depending upon which node the list is representing.
#' @export
#'
build_single_tree <- function(X, Y, mtree, min_leaf,V_inv,Command){
  model=rep( list(NULL), 10000 )
  i=1
  Index=1:nrow(X)
  
  model=split_node(X,Y,mtree,Index,i,model,min_leaf,V_inv,Command)
  return(model)
}