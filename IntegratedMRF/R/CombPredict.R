#' Integrated Prediction of Testing samples using Combination Weights from integrated RF or MRF model 
#' 
#' Generates Random Forest or Multivariate Random Forest model for each subtype of dataset and predicts testing samples using the generated models. 
#' Subsequently, the prediction for different subtypes of dataset are combined using the Combination weights generated from 'Combination' function.
#'  
#' @param finalX List of Matrices where each matrix represent a specific data subtype (such as genomic characterizations for 
#' drug sensitivity prediction). Each subtype can have different types of features. For example, if there are three subtypes containing
#'  100, 200 and 250 features respectively,  finalX will be a list containing 3 matrices of sizes M x 100, M x 200 and M x 250 
#'  where M is the number of Samples.
#' @param finalY_train A M x T matrix of output features for training samples, where M is number of samples and
#' T is the number of output features. The dataset is assumed to contain no missing values. If there are missing values, an imputation method 
#' should be applied before using the function. A function 'Imputation' is included within the package.
#' @param Cell It contains a list of samples (the samples can be represented either numerically by indices or by names) for each data subtype. 
#' For the example of 3 data subtypes, it will be a list containing 3 arrays where each array contains the sample information for each data subtype.
#' @param finalY_train_cell Cell lines of output features for training samples
#' @param finalY_test_cell Cell lines of output features for testing samples
#' @param n_tree number of trees in the forest
#' @param m_feature Number of randomly selected features considered for a split in each regression tree node.
#' @param min_leaf minimum number of samples in the leaf node 
#' @param Serial Consists of a  list of all combinations of different subtypes of a dataset (except for the case with no dataset being selected). 
#' For example, if a 
#' dataset has 3 subtypes, then Serial is a list of size 2^3-1=7.  The ordering of the seven sets will be [1 2 3], [1 2], [1 3], [2 3], [1], [2], [3]
#' @param Coeff Combination Weights. The user can supply the weights based on either Bootstrap, Re-substitution, 0.632Bootstrap or Leave-one-out 
#' error estimation approaches.
#' 
#' @return Final Prediction of testing samples based on provided testing sample names.
#' @examples
#' ## This example is the extension of the example given in the "Combination" function.
#' ##finalX, finalY_train, Cell, finalY_train_cell, n_tree,m_feature, 
#' ##min_leaf, Serial are the same variable provided for the "Combination" function.
#' ##The Weight is one of the combination weight(Bootstrap, Re-substituition, 
#' ##0.632Bootstrap or Leave-one-out) got from the result of "Combination" function.
#' #Uncomment the following line to run the code
#' #Drug_Sen_test <- read.xlsx("Output_Response_Test.xlsx", sheet = 1,  startRow = 1, colNames = TRUE)
#' #finalY_test_cell=Drug_Sen_test[,1]
#' #Final_Prediction=CombPredict(finalX,finalY_train,Cell,
#' #+finalY_train_cell,finalY_test_cell,n_tree,m_feature,min_leaf,Serial,Weight)
#' @details
#' Input matrix and output response of training samples have been used to build Random Forest or Multivariate Random Forest model for each subtype of 
#' a dataset. These models are used to calculate prediction of 
#' testing samples for each subtype separately. Subsequently Combination Weights (different errors have different combination weights
#'  and the user should select the one to be used) are used to integrate the predictions from data subtypes. 
#'  Note that the combination weights are linear regression coefficients generated using the training samples. 
#'  
#'  The specific set of combination weights to be used for testing samples will depend on the number of data subtypes available 
#'  for the testing samples. Note that not all subtype information maybe available for all samples. 
#'  As an example with three data subtypes, a testing sample with all subtype data available will using 
#'  the combination weights corresponding to Serial [1 2 3] where if subtype 3 is not available, the function will 
#'  using the combination weights corresponding to Serial [1 2].
#' @export
CombPredict <- function(finalX,finalY_train,Cell,finalY_train_cell,finalY_test_cell,n_tree,m_feature,min_leaf,Serial,Coeff){
  Common_cell_dataset=NULL
  for (q in 1:length(Cell)){
    Common_cell_dataset[[q]]=intersect(finalY_train_cell,Cell[[q]])
  }
  Variable_number=ncol(finalY_train)
  if (Variable_number>1){
    Command=2
  }else if(Variable_number==1){
    Command=1
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
  
  Y_hat_test=NULL
  final_test=NULL
  Common_cell_test=NULL
  for (q in 1:length(Cell)){
    Common_cell_test[[q]]=intersect(finalY_test_cell,Cell[[q]])
    Cell_ind=match(Common_cell_test[[q]],Cell[[q]])
    final_test[[q]]=finalX[[q]][Cell_ind, ]
    final_test[[q]]=matrix(as.numeric(final_test[[q]]),nrow = dim(final_test[[q]])[1], ncol = dim(final_test[[q]])[2])
    
    Y_hat_test[[q]] = build_forest_predict(final_dataset[[q]],finalY_dataset[[q]], n_tree, m_feature, min_leaf, final_test[[q]])
  }
  Common_cell=NULL
  Final_test=NULL
  Drug_sensitivity_cell_test2=finalY_test_cell
  for (q in 1:length(Serial)){
    D=NULL
    D=Drug_sensitivity_cell_test2
    for (check in 1:length(Serial[[q]])){
      D=intersect(D,Cell[Serial[[q]][check]][[1]])
    }
    Common_cell[[q]]=D
    Drug_sensitivity_cell_test2=setdiff(Drug_sensitivity_cell_test2,D)
    
    F_test=matrix(rep(0,length(Common_cell[[q]])*Variable_number),ncol=Variable_number)
    match_ind=NULL
    W=Serial[[q]]
    for (RR in 1:Variable_number){
      for (R in 1:length(W)){
        match_ind[[R]]=match(Common_cell[[q]],Common_cell_test[[W[R]]])
        F_test[,RR]=F_test[,RR]+Coeff[[q]][R,RR]*matrix(Y_hat_test[[W[R]]][match_ind[[R]],RR],ncol=1)
      }
    }
    
    Final_test[[q]]=F_test
  }
  Final_result=NULL
  for (q in 1:length(Serial)){
    if (length(Common_cell[[q]])>0){
      Final_result=rbind(Final_result,cbind(matrix(Common_cell[[q]],ncol=1),matrix(Final_test[[q]],ncol=Variable_number)))
    }
  }
  match_ind=match(finalY_test_cell,Final_result[,1])
  Final_prediction=Final_result[match_ind,]
  return(Final_prediction)
}