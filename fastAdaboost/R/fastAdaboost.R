#' fastAdaboost: fast adaboost implementation for R
#' 
#' fastAdaboost provides a blazingly fast implementation of both discrete 
#' and real adaboost algorithms, based on a C++ backend. The goal of the 
#' package is to provide fast performance for large in-memory data sets.
#' 
#'@examples
#'fakedata <- data.frame( X=c(rnorm(100,0,1),rnorm(100,1,1)), Y=c(rep(0,100),rep(1,100) ) )
#'fakedata$Y <- factor(fakedata$Y)
#'test_adaboost <- adaboost(Y~X, fakedata, 10)
#'pred <- predict( test_adaboost,newdata=fakedata)
#'print(pred$error)
#'
#'@examples
#'fakedata <- data.frame( X=c(rnorm(100,0,1),rnorm(100,1,1)), Y=c(rep(0,100),rep(1,100) ) )
#'fakedata$Y <- factor(fakedata$Y)
#'test_real_adaboost <- real_adaboost(Y~X, fakedata, 10)
#'pred <- predict(test_real_adaboost,newdata=fakedata)
#'print(pred$error)
#'
#'@references 
#'Freund, Y. and Schapire, R.E. (1996):\dQuote{Experiments with a new boosting algorithm}
#'. \emph{In Proceedings of the Thirteenth International Conference on Machine Learning}, 
#'pp. 148--156, Morgan Kaufmann.
#'
#'Zhu, Ji, et al. \dQuote{Multi-class adaboost} \emph{Ann Arbor} 1001.48109 (2006): 1612.
#' @docType package
#' @author Sourav Chatterjee
#' @name fastAdaboost
#' 
#' @useDynLib fastAdaboost
#' @importFrom Rcpp sourceCpp
NULL