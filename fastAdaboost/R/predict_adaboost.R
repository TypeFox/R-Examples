#'predict method for adaboost objects
#'
#'predictions for model corresponding to adaboost.m1 algorithm
#'
#'makes predictions for an adaboost object on a new dataset.
#'The target variable is not required 
#'for the prediction to work.
#'However, the user must ensure that the test data has the same 
#'columns which were used as inputs to fit the original model.
#'The error component of the prediction object(as in 
#'\code{pred$error}) can be used to get the error of the 
#'test set if the test data is labeled.
#'
#'@seealso \code{\link{adaboost}}

#'@import rpart
#'@param object an object of class adaboost
#'@param newdata dataframe on which we are looking to predict
#'@param ... arguments passed to predict.default
#'@return predicted object, which is a list with the following components
#'\item{formula}{the formula used.}
#'\item{votes}{total weighted votes achieved by each class}
#'\item{class}{the class predicted by the classifier}
#'\item{prob}{a matrix with predicted probability of each class for each observation}
#'\item{error}{The error on the test data if labeled, otherwise \code{NA}}
#'@export
#'@examples
#'fakedata <- data.frame( X=c(rnorm(100,0,1),rnorm(100,1,1)), Y=c(rep(0,100),rep(1,100) ) )
#'fakedata$Y <- factor(fakedata$Y)
#'test_adaboost <- adaboost(Y~X, fakedata, 10)
#'pred <- predict( test_adaboost,newdata=fakedata)
#'print(pred$error)
#'print( table(pred$class,fakedata$Y) )

predict.adaboost <- function(object, newdata,...)
{
  return( predict_adaboost_internal(object,newdata,"M1") )
}



#'predict method for real_adaboost objects
#'
#'predictions for model corresponding to real_adaboost algorithm
#'
#'makes predictions for an adaboost object on a new dataset
#'using the real_adaboost algorithm.
#' The target variable is not required 
#'for the prediction to work.
#'However, the user must ensure that the test data has the same 
#'columns which were used as inputs to fit the original model.
#'The error component of the prediction object(as in 
#'\code{pred$error}) can be used to get the error of the 
#'test set if the test data is labeled.
#'
#'@import rpart
#'@param object an object of class real_adaboost
#'@param newdata dataframe on which we are looking to predict
#'@param ... arguments passed to predict.default
#'@return predicted object, which is a list with the following components
#'\item{formula}{the formula used.}
#'\item{votes}{total weighted votes achieved by each class}
#'\item{class}{the class predicted by the classifier}
#'\item{prob}{a matrix with predicted probability of each class for each observation}
#'\item{error}{The error on the test data if labeled, otherwise \code{NA}}
#'@export
#'@examples
#'fakedata <- data.frame( X=c(rnorm(100,0,1),rnorm(100,1,1)), Y=c(rep(0,100),rep(1,100) ) )
#'fakedata$Y <- factor(fakedata$Y)
#'test_real_adaboost <- real_adaboost(Y~X, fakedata, 10)
#'pred <- predict(test_real_adaboost,newdata=fakedata)
#'print(pred$error)
#'print( table(pred$class,fakedata$Y) )
#'
#'@seealso \code{\link{real_adaboost}}
#'
predict.real_adaboost <- function(object, newdata,...)
{
  return( predict_adaboost_internal(object,newdata,"SAMME.R") ) 
}





#'internal method which predicts both adaboost_M1 and 
#'real adaboost
#'@noRd
#'@param object an object of class adaboost
#'@param newdata dataframe on which we are looking to predict
#'@param method M1 or SAMME.R
#'@keywords internal
predict_adaboost_internal <-function(object, newdata, method)
{
  if(!all(method %in% c("M1","SAMME.R") ))
    stop(paste("method must be M1 or SAMME.R. It is",method), call.=F)
  formula <- object$formula
  tree_list <- object$trees
  coeff_vector <- object$weights
  classnames_map <- object$classnames
  num_examples <- nrow(newdata)
  if(method == "M1")
  {
    cpp_list <- predict_adaboost_(tree_list, coeff_vector, newdata,
                                  num_examples, wrap_rpart_predict,classnames_map)
  }
  else
  {
    cpp_list <- predict_real_adaboost_(tree_list, coeff_vector, newdata,
                                       num_examples, wrap_rpart_predict_real)
  }
  votes <- cpp_list$votes
  predicted_class_int <- cpp_list$class #this is 0 or 1
  predicted_class <- ifelse(predicted_class_int == 0, classnames_map["A"],classnames_map["B"])
  predicted_class <- factor(predicted_class)
  prob_mat = cpp_list$prob
  
  #if data is labeled, calculate prediction error
  test_error <- NA
  depvar_name <- object$dependent_variable 
  if( depvar_name %in% names(newdata))
  {
    vardep = ifelse(newdata[,depvar_name]==classnames_map["A"],0,1)
    test_error <- calculate_test_error_(vardep, predicted_class_int)
  }
    
  predictor <- list(formula = formula, votes = votes, 
                    class = predicted_class,prob = prob_mat,
                    error = test_error)
  
  
  return(predictor)
}

