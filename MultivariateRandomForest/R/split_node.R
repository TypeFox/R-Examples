#' Splitting Criteria of all the nodes of the tree
#' 
#' Stores the Splitting criteria of all the nodes of a tree in a list 
#'  
#' @param X Input Training matrix of M x N, M is the number of training samples and N is the number of features
#' @param Y Output Training response of M x T, M is the number of samples and T is number of ouput Features(Response)
#' @param mtree Number of randomly selected features used for each split
#' @param Index Index of training samples
#' @param i Number of split. Used as an index, which indicates where in the list the splitting
#' criteria of this split will be stored.
#' @param model A list of lists with the spliting criteria of all the split of the nodes. In each iteration,
#' a new list is included with the spliting criteria of the new split of a node.
#' @param min_leaf Minimum number of samples in the leaf node. If a node has less than equal to min_leaf samples,
#'  then there will be no splitting in that node and this node is a leaf node.
#' @param V_inv Covariance matrix of Output Feature matrix
#' @param Command 1 for RF and 2 for MRF depending on the method
#' @return Model A list of lists with the spliting criteria of all the split of the nodes. In each iteration,
#' the Model is updated with a new list which includes the spliting criteria of the new split of a node.
#' @details
#' This function calculates the spltting criteria of a node and stores the information in a list format. If the 
#' node is a parent node, then index of left and right nodes and feature number and threshold value of the feature
#' for the split has been stored. While if the node is a leaf, then output feature matrix of the samples for
#' the node has been stored as a list in the Model.
#' @export
split_node <- function(X,Y,mtree,Index,i,model,min_leaf,V_inv,Command){
  ii=NULL
  Index_left=NULL
  Index_right=NULL
  if(length(Index)>min_leaf){ #create problem with 2
    Result = spliting(X,Y,mtree,Index,V_inv,Command)
    Index_left=Result[[1]]
    Index_right=Result[[2]]
    if(i==1){
      Result[[5]]=c(2,3)
    }else{
      j=1
      while (length(model[[j]])!=0){
        j=j+1
      }
      Result[[5]]=c(j,j+1)
    }
    
    model[[i]]=Result
    k=i
    i=1 #maybe unnecessary
    while (length(model[[i]])!=0){
      i=i+1
    } 
    model[[Result[[5]][1]]]=Result[[1]]
    model[[Result[[5]][2]]]=Result[[2]]
    
    model=split_node(X,Y,mtree,Index_left,model[[k]][[5]][1],model,min_leaf,V_inv,Command)
    
    model=split_node(X,Y,mtree,Index_right,model[[k]][[5]][2],model,min_leaf,V_inv,Command)
    
    
  }else{
    ii[[1]]=matrix(Y[Index,],ncol=ncol(Y))
    model[[i]]=ii
  }
  
  
  return(model)
}