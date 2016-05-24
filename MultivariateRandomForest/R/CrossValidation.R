#' Matrix of Input and Output of Cross validation
#' 
#' Generates Cross validated Input Matices and Output Vectors, where number of fold in cross validation is user defined
#'  
#' @param X M x N Input matrix, M is the number of samples and N is the number of features
#' @param Y output Response as column vector 
#' @param F Number of Fold in cross validation
#' @return List with the following components: 
#' \item{TrainingData}{List of matices with matrix containing Cross Validates Training Data, where number of list
#' equal to user defined cross validation}
#' \item{TestingData}{List of matices with matrix containing Cross Validates Testing Data, where number of list
#' equal to user defined cross validation}
#' \item{OutputTrain}{List of matices with matrix containing Cross Validates Training Response Data, where number of list
#' equal to user defined cross validation}
#' \item{OutputTest}{List of matices with matrix containing Cross Validates Testing Response Data, where number of list
#' equal to user defined cross validation}
#' \item{FoldedIndex}{Index of Different Fold}
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