#' @title OASIS Training
#' @description This function trains the OASIS model from a \code{data.frame} 
#' produced by the function \code{\link{oasis_train_dataframe}} 
#' @param ... \code{data.frame}(s) produced by the 
#' \code{\link{oasis_train_dataframe}} function 
#' @param remove_preproc a logical stating if preprocessed volumes 
#' need to be removed from the \code{data.frame}s 
#' @import fslr
#' @export 
#' @return Returns a \code{glm} object containing the trained OASIS 
#' coefficients to be used by the function \code{\link{oasis_predict}}.
#' @examples \dontrun{
#' my_oasis_model <- oasis_training(oasis_dataframe_1, oasis_dataframe_2)  
#' }
#' @import oro.nifti
oasis_training <- function(..., ##dataframes from function 
                           remove_preproc  = FALSE) 
{
  list_of_train_dataframes <- list(...)
  if (remove_preproc  == TRUE) {
    list_of_train_dataframes  <- lapply( list_of_train_dataframes, function(x) x[[1]])
  }
  train_vectors_multi <- do.call(rbind, list_of_train_dataframes)  
  train_vectors_multi <- as.data.frame(train_vectors_multi)
  ##fit the oasis model 
  oasis_model <- glm(formula = GoldStandard ~ FLAIR_10 *FLAIR  +
                      FLAIR_20*FLAIR + PD_10 *PD  + PD_20 *PD +
                      T2_10 *T2 +  T2_20 *T2 + T1_10 *T1 +
                      T1_20 *T1, 
                     data = train_vectors_multi, 
                     family = binomial)
  
  ##clean up the oasis model 
  oasis_model$y = c()
  oasis_model$model = c()
  oasis_model$residuals = c()
  oasis_model$fitted.values = c()
  oasis_model$effects = c()
  oasis_model$qr$qr = c()  
  oasis_model$linear.predictors = c()
  oasis_model$weights = c()
  oasis_model$prior.weights = c()
  oasis_model$data = c()
  oasis_model$family$variance = c()
  oasis_model$family$dev.resids = c()
  oasis_model$family$aic = c()
  oasis_model$family$validmu = c()
  oasis_model$family$simulate = c()
  attr(oasis_model$terms,".Environment") = c()
  attr(oasis_model$formula,".Environment") = c()
  
  return(oasis_model)
}