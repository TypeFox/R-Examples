#' Prediction of testing sample in a node
#' 
#' Predicts a testing sample in a node, means which child node it will go using the splitting criteria of the model
#' or, prediction results if the node is leaf. 
#'  
#' @param Single_Model Model of a particular tree
#' @param i Number of split. Used as an index, which indicates where in the list the splitting
#' criteria of this split has been stored.
#' @param X_test Testing samples of Q x N, Q is the number of testing samples and N is the number of features(same order and
#' size used as training) 
#' @param Variable_number Number of Output Features
#' @return Prediction result of a Testing samples in a node
#' @details
#' In a particular node, a testing sample can go in two directions, left node and right node. If the value of splitting feature
#' of the testing sample is less than thresohold splitting value, it will go to the left otherwise in the right node.
#' While if the node is leaf, then average of output features are returned as the prediction results.
#' @export
predicting <- function(Single_Model,i,X_test,Variable_number){
  
  Result=NULL
  
  if(length(Single_Model[[i]])==5){
    feature_no=Single_Model[[i]][[3]]
    feature_value=X_test[feature_no]
    if(feature_value<Single_Model[[i]][[4]]){  #feature value less than threshold value
      #i=i*2+1
      Result=predicting(Single_Model,Single_Model[[i]][[5]][1],X_test,Variable_number)
    }else{                                    #feature value greater than threshold value
      #i=i*2+2
      Result=predicting(Single_Model,Single_Model[[i]][[5]][2],X_test,Variable_number)
    }
  }else{
    Result=matrix(  0*Variable_number,  ncol=Variable_number)
    if (Variable_number>1){
      for (jj in 1:Variable_number){
        Result[,jj]=mean(Single_Model[[i]][[1]][,jj])
      }
    }else {
      for (jj in 1:Variable_number){
        Result[,jj]=mean(unlist(Single_Model[[i]][[1]]))
      }
    }
    
  }
  return(Result)
}