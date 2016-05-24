#'Adaboost.M1 algorithm
#'
#'Implements Freund and Schapire's Adaboost.M1 algorithm
#'
#'This implements the Adaboost.M1 algorithm for a binary classification task.
#'The target variable must be a factor with exactly two levels.
#'The final classifier is a linear combination of weak decision tree classifiers.
#'
#'@import rpart
#'@param formula Formula for models
#'@param data Input dataframe
#'@param nIter no. of classifiers 
#'@param ... other optional arguments, not implemented now
#'@return object of class adaboost
#'@export 
#'@examples 
#'fakedata <- data.frame( X=c(rnorm(100,0,1),rnorm(100,1,1)), Y=c(rep(0,100),rep(1,100) ) )
#'fakedata$Y <- factor(fakedata$Y)
#'test_adaboost <- adaboost(Y~X, data=fakedata,10)
#'@seealso \code{\link{real_adaboost}}, \code{\link{predict.adaboost}}
#'@references Freund, Y. and Schapire, R.E. (1996):\dQuote{Experiments with a new boosting algorithm}
#'. \emph{In Proceedings of the Thirteenth International Conference on Machine Learning}, 
#'pp. 148--156, Morgan Kaufmann.

adaboost <-function(formula, data, nIter,...)
{
  theCall <- match.call()
  adaboost_object <- adaboost_fast(formula,data,nIter, method="M1")
  adaboost_object$call <- theCall
  return(adaboost_object)
}

