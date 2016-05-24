#'Real Adaboost algorithm
#'
#'Implements Zhu et al's real adaboost or SAMME.R algorithm
#'
#'This implements the real adaboost algorithm for a binary classification 
#'task. The target variable must be a factor with exactly two levels.
#'The final classifier is a linear combination of weak decision tree
#'classifiers. Real adaboost uses the class probabilities of the weak classifiers
#'to iteratively update example weights. It has been found to have lower
#'generalization errors than adaboost.m1 for the same number of iterations.
#'
#'@import rpart
#'@param formula Formula for models
#'@param data Input dataframe
#'@param nIter no. of classifiers 
#'@param ... other optional arguments, not implemented now
#'@return object of class real_adaboost
#'@export 
#'@examples 
#'fakedata <- data.frame( X=c(rnorm(100,0,1),rnorm(100,1,1)), Y=c(rep(0,100),rep(1,100) ) )
#'fakedata$Y <- factor(fakedata$Y)
#'test_adaboost <- real_adaboost(Y~X, data=fakedata,10)
#'@seealso \code{\link{adaboost}},\code{\link{predict.real_adaboost}}
#'@references Zhu, Ji, et al. \dQuote{Multi-class adaboost} \emph{Ann Arbor} 1001.48109 (2006): 1612.


real_adaboost <-function(formula, data, nIter,...)
{
  theCall <- match.call()
  adaboost_object <- adaboost_fast(formula,data,nIter, method="SAMME.R")
  adaboost_object$call <- theCall
  return(adaboost_object)
}

