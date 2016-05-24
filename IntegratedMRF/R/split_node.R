#' Splitting Criteria of all the nodes of the tree
#' 
#' Stores the Splitting criteria of all the nodes of a tree in a list 
#'  
#' @param X Input Training matrix of size M x N, M is the number of training samples and N is the number of features
#' @param Y Output Training response of size M x T, M is the number of samples and T is the number of output responses
#' @param m_feature Number of randomly selected features considered for a split in each regression tree node
#' @param Index Index of training samples
#' @param i Number of split. Used as an index, which indicates where in the list the splitting
#' criteria of this split will be stored.
#' @param model A list of lists with the spliting criteria of all the node splits. In each iteration,
#' a new list is included with the spliting criteria of the new split of a node.
#' @param min_leaf Minimum number of samples in the leaf node. If a node has less than or, equal to min_leaf samples,
#'  then there will be no splitting in that node and the node is a leaf node.
#' @param Cov_Y Covariance matrix of Output Response matrix for MRF(Give Zero for RF)
#' @param Command 1 for univariate Regression Tree (corresponding to RF) and 2 for Multivariate Regression Tree (corresponding to MRF)
#' @return Model: A list of lists with the splitting criteria of all the split of the nodes. In each iteration, the Model is 
#' updated with a new list that includes the splitting criteria of the new split of a node.
#' @details
#' This function calculates the splitting criteria of a node and stores the information in a list format. 
#' If the node is a parent node, then indices of left and right nodes and feature number and threshold value 
#' of the feature for the split are stored. If the node is a leaf, the output feature matrix of the samples 
#' for the node are stored as a list.
#' @export
split_node <- function(X,Y,m_feature,Index,i,model,min_leaf,Cov_Y,Command){
  ii=NULL
  Index_left=NULL
  Index_right=NULL
  if(length(Index)>min_leaf){ #create problem with 2
    Result = spliting(X,Y,m_feature,Index,Cov_Y,Command)
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
    
    model=split_node(X,Y,m_feature,Index_left,model[[k]][[5]][1],model,min_leaf,Cov_Y,Command)
    
    model=split_node(X,Y,m_feature,Index_right,model[[k]][[5]][2],model,min_leaf,Cov_Y,Command)
    
    
  }else{
    ii[[1]]=matrix(Y[Index,],ncol=ncol(Y))
    model[[i]]=ii
  }
  
  
  return(model)
}