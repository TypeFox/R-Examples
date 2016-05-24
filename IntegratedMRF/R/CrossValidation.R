#' Generate training and testing samples for cross validation
#' 
#' Generates Cross Validation Input Matrices and Output Vectors for training and testing, where number of folds in cross validation is user defined.
#'  
#' @param X M x N Input matrix, M is the number of samples and N is the number of features
#' @param Y output Response as column vector 
#' @param F Number of Folds
#' @return List with the following components: 
#' \item{TrainingData}{List of matrices where each matrix contains a fold of Cross Validation Training Data, 
#' where the number of matrices is equal to F}
#' \item{TestingData}{List of matrices where each matrix contains a fold of Cross Validation Testing Data, 
#' where the number of matrices is equal to F}
#' \item{OutputTrain}{List of matrices where each matrix contains a fold of Cross Validation Training Output Feature Data, 
#' where the number of matrices is equal to F}
#' \item{OutputTest}{List of matrices where each matrix contains a fold of Cross Validation Testing Output Feature Data, 
#' where the number of matrices is equal to F}
#' \item{FoldedIndex}{Index of Different Fold. (e.g., for Sample Index 1:6 and 3 fold, FoldedIndex are [1 2 3 4], [1 2 5 6], [3 4 5 6])}
#' @export

CrossValidation <- function(X,Y,F) {
  
  DrugNumber=nrow(Y)
  Index=1:DrugNumber
  
  FoldedIndex=NULL
  FF=F-1
  for (i in 1:FF) {
    FoldedIndex[[i]]=sample(Index,floor(DrugNumber/F))
    Index=setdiff(Index,FoldedIndex[[i]])
    if (i==FF) {
      FoldedIndex[[i+1]]=Index
    }
  }
  TestingIndex=NULL
  TrainingIndex=NULL
  TrainingData=NULL
  TestingData=NULL
  OutputTrain=NULL
  OutputTest=NULL
  for (i in 1:F){
    TestingIndex[[i]]=FoldedIndex[[i]]
    TrainingIndex[[i]]=setdiff(1:DrugNumber,TestingIndex[[i]])
    ## Input 
    TrainingData[[i]]=X[TrainingIndex[[i]],,drop=FALSE]  #rows of training data of matrix X
    TestingData[[i]]=X[TestingIndex[[i]],,drop=FALSE]      #rows of testing data of matrix X
    ## Output
    OutputTrain[[i]]=Y[TrainingIndex[[i]],,drop=FALSE]        #rows of training data of matrix Y
    OutputTest[[i]]=Y[TestingIndex[[i]],,drop=FALSE]           #rows of testing data of matrix Y
  }
  result=list(TrainingData,TestingData,OutputTrain,OutputTest,FoldedIndex)
  return(result)
}