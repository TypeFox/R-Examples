#' Prediction for testing samples using specific combination weights from integrated RF or MRF model 
#' 
#' Generates Random Forest (One Output Feature) or Multivariate Random Forest (More than One Output Feature) 
#' model for each subtype of dataset and predicts testing samples using these models. The predictions are 
#' combined using the specific combination weights provided by the user. For the input combination weights, 
#' the testing cell lines should have the subtype data corresponding to the non-zero weight subtypes. 
#'  
#' @param finalX List of Matrices where each matrix represent a specific data subtype (such as genomic characterizations for 
#' drug sensitivity prediction). Each subtype can have different types of features. For example, if there are three subtypes containing
#'  100, 200 and 250 features respectively,  finalX will be a list containing 3 matrices of sizes M x 100, M x 200 and M x 250 
#'  where M is the number of Samples.
#' @param finalY_train A M x T matrix of output features for training samples, where M is number of samples and T is the number of output features. 
#' The dataset is assumed to contain no missing values. If there are missing values, an imputation method should be applied before using the function. 
#' A function 'Imputation' is included within the package.
#' @param Cell It contains a list of samples (the samples can be represented either numerically by indices or by names) for each data subtype. 
#' For the example of 3 data subtypes, it will be a list containing 3 arrays where each array contains the sample information for each data subtype.
#' @param finalY_train_cell Sample names of output features for training samples
#' @param finalY_test_cell Sample names of output features for testing samples(All these testing samples
#' must have features for each subtypes of dataset)
#' @param n_tree Number of trees in the forest
#' @param m_feature Number of randomly selected features considered for a split in each regression tree node
#' @param min_leaf Minimum number of samples in the leaf node 
#' @param Serial Consists of a  list of all combinations of different subtypes of a dataset (except for the case with no dataset being selected). 
#' For example, if a 
#' dataset has 3 subtypes, then Serial is a list of size 2^3-1=7.  The ordering of the seven sets will be [1 2 3], [1 2], [1 3], [2 3], [1], [2], [3]
#' @param Coeff Combination Weights (user defined or some combination weights generated using the 'Combination' function). 
#' The size must be C, which is equal to the number of subtypes of dataset given in finalX.
#' 
#' @return Final Prediction of testing samples based on provided testing sample names
#' @examples
#' ## This example is the extension of the example given in the "Combination" function.
#' ##finalX, finalY_train, Cell, finalY_train_cell, n_tree,m_feature, 
#' ##min_leaf, Serial are the same variable provided for the "Combination" function.
#' ##The Weight is a vector of size length(Cell)x1, which is user defined.
#' #Uncomment the following line to run the code
#' #Drug_Sen_test <- read.xlsx("Output_Response_Test.xlsx", sheet = 1,  startRow = 1, colNames = TRUE)
#' #finalY_test_cell=Drug_Sen_test[,1]
#' #Final_Prediction=CombPredictSpecific(finalX,finalY_train,Cell,
#' #+finalY_train_cell,finalY_test_cell,n_tree,m_feature,min_leaf,Serial,runif(length(Cell)*1))
#' @details
#' Input feature matrix and output feature matrix have been used to generate Random Forest (One Output Feature) 
#' or Multivariate Random Forest (More than One Output Feature) model for each subtype of dataset separately. 
#' The prediction of testing samples using each subtype trained model is generated. The predictions are combined 
#' using the specific combination weights provided by the user. For the input combination weights, the testing cell lines 
#' should have the subtype data corresponding to the non-zero weight subtypes. For instance, if combination weights is 
#' [0.6 0.3 0 0.1], then the subtype 1, 2 and 4 needs to be present for the testing samples. Furthermore, all the features 
#' should be present for the required subtypes for the testing samples.
#' @export
CombPredictSpecific <- function(finalX,finalY_train,Cell,finalY_train_cell,finalY_test_cell,n_tree,m_feature,min_leaf,Serial,Coeff){
  Common_cell_dataset=NULL
  for (q in 1:length(Cell)){
    Common_cell_dataset[[q]]=intersect(finalY_train_cell,Cell[[q]])
  }
  Variable_number=ncol(finalY_train)
  if (ncol(finalY_train)==1){
    Command=1
  }else if (ncol(finalY_train)>1){
    Command=2
  }
  
  final_dataset=NULL
  finalY_dataset=NULL
  for (q in 1:length(Cell)){
    Cell_ind=match(Common_cell_dataset[[q]],Cell[[q]])
    final_dataset[[q]]=finalX[[q]][Cell_ind, ]
    final_dataset[[q]]=matrix(as.numeric(final_dataset[[q]]),nrow = dim(final_dataset[[q]])[1], ncol = dim(final_dataset[[q]])[2])
    
    Cell_ind_Y=match(Common_cell_dataset[[q]],finalY_train_cell)
    finalY_dataset[[q]]=matrix(finalY_train[Cell_ind_Y,],ncol=Variable_number)
  }
  
  Y_hat_test=NULL#matrix(rep(0,length(Cell)*length(finalY_test_cell)),nrow=length(finalY_test_cell))
  final_test=NULL
  Common_cell_test=NULL
  for (q in 1:length(Cell)){
    Common_cell_test=intersect(finalY_test_cell,Cell[[q]])
    Cell_ind=match(Common_cell_test,Cell[[q]])
    final_test[[q]]=finalX[[q]][Cell_ind, ]
    final_test[[q]]=matrix(as.numeric(final_test[[q]]),nrow = dim(final_test[[q]])[1], ncol = dim(final_test[[q]])[2])
    
    Y_hat_test[[q]] = build_forest_predict(final_dataset[[q]],finalY_dataset[[q]], n_tree, m_feature, min_leaf, final_test[[q]])
  }
  Y_hat_temp=matrix(rep(0,length(Cell)*length(finalY_test_cell)),nrow=length(finalY_test_cell))
  Final_prediction=matrix(Common_cell_test,ncol=1)
  for (RR in 1:Variable_number){
    for (q in 1:length(Cell)){
      Y_hat_temp[,q]=Y_hat_test[[q]][,RR]
    }
    Final_prediction=cbind(Final_prediction,Y_hat_temp%*%matrix(Coeff,nrow=length(Cell)))
  }
  return(Final_prediction)
}