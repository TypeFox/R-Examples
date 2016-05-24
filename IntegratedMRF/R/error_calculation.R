#' Error calculation for integrated model 
#' 
#' Combines Prediction from different data subtypes through Least Square Regression and computes Mean Absolute Error, 
#' Mean Square Error and Pearson Correlation Coefficient between Integrated Prediction and Original Output feature.
#'  
#' @param final_pred A n x p matrix of predicted features, where n is the number of samples and p is the number of data subtypes with prediction
#' @param final_actual A n x 1 vector of original output responses
#' @return List with the following components: 
#' \item{Integrated Prediction}{Integrated Prediction based on combining predictions from data subtypes using Least Square Regression}
#' \item{error_mae}{Mean Absolute Error between Integrated Prediction and Original Output Responses}
#' \item{error_mse}{Mean Square Error between Integrated Prediction and Original Output Responses}
#' \item{error_corr}{Pearson Correlation Coefficient between Integrated Prediction and Original Output Responses}
#' @details 
#' If final_pred is a vector, it refers to the prediction result for one subtype of dataset and this function will return 
#' Mean Absolute Error, Mean Square Error and Pearson Correlation Coefficient between between predicted and Original Output response. 
#' If final_pred is a matrix containing prediction results for more than one subtype of dataset, Least Square 
#' Regression will be using to calculate the weights to combine the predictions and generate an integrated prediction of size n x 1. 
#' Subsequently, Mean Absolute Error, Mean Square Error and Pearson Correlation Coefficient between 
#' Integrated Prediction and Original Output response are calculated.
#' @seealso \code{lsei}
#' @export
error_calculation <- function(final_pred,final_actual){
  #   library(limSolve)
  a0=matrix(limSolve::lsei(A=final_pred, B=final_actual, E=rep(1,dim(final_pred)[2]), F=1)$X, ncol=1)
  final_pred2=final_pred%*%a0
  
  k=0
  err1=NULL
  err2=NULL
  for (j in 1:length(final_actual)){
    if(final_actual[j,]>=1e-3){
      k=k+1
      err1[k]=abs(final_pred2[j,]-final_actual[j,])
      err2[k]=(final_pred2[j,]-final_actual[j,])^2
      #error1(Drug)=sum(abs((final_gen_char-final_sensi)./final_sensi))/length(final_sensi);
    }
  }
  error=NULL
  error_mae=sum(err1)/k
  error_mse=sum(err2)/k
  corr=stats::cor(final_actual,final_pred2)  
  error[[1]]=final_pred2
  error[[2]]=error_mae
  error[[3]]=error_mse
  error[[4]]=corr
  return(error)
}