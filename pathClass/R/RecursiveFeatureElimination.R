#' Recursive Feature Elimination (RFE)
#'
#' Implementation of the Recursive Feature Elimination (RFE) algorithm.
#'
#' @param x a p x n matrix of expression measurements with p samples and n genes.
#' @param y a factor of length p comprising the class labels.
#' @param DEBUG should debugging information be plotted.
#' @param scale a character vector defining if the data should be centered and/or scaled.
#' Possible values are \emph{center} and/or \emph{scale}. Defaults to \code{c('center', 'scale')}.
#' @param Cs soft-margin tuning parameter of the SVM. Defaults to \code{10^c(-3:3)}.
#' @param stepsize amount of features that are discarded in each step of the feature elimination. Defaults to 10\%.
#' @return a RFE fit object.
#' \code{features} = selected features
#' \code{error.bound} = span bound of the model
#' \code{fit} = fitted SVM model
#' @export
#' @callGraphPrimitives
#' @note The optimal number of features is found by using the span estimate. See Chapelle, O., Vapnik, V., Bousquet, O., and Mukherjee, S. (2002). Choosing multiple parameters for support vector machines. \emph{Machine Learning}, 46(1), 131-159.
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' \dontrun{
#' library(Biobase)
#' data(sample.ExpressionSet)
#' x <- t(exprs(sample.ExpressionSet))
#' y <- factor(pData(sample.ExpressionSet)$sex)
#' res.rfe <- crossval(x,y,DEBUG=TRUE,theta.fit=fit.rfe,folds=2,repeats=1,parallel=TRUE,Cs=10^(-3:3))
#' }
fit.rfe = function(x, y, DEBUG=FALSE, scale=c('center', 'scale'), Cs=10^c(-3:3), stepsize=0.1){

  best.bound = Inf
  feat = colnames(x)

  while(NCOL(x) > 1){		
    if(DEBUG) cat(NCOL(x),' Features left.\n')
    fit = svm.fit(x=x, y=y, Cs=Cs, scale=scale, DEBUG=DEBUG)
    
    if(fit$error.bound <= best.bound){
      best.bound = fit$error.bound
      feat  = colnames(x)			
      best.fit   = fit
      if(DEBUG) cat('Model Updated. Spanbound=',best.bound,', C=',best.fit$C,', ', length(feat),'features.\n')
    }
    
    ord      = order(fit$w)
    remove   = colnames(x)[ord[1:round(NCOL(x)*stepsize)]]	# make it quicker!
    x        = x[, setdiff(colnames(x), remove), drop=FALSE]						
  }
  if(DEBUG) cat('Best Model is: Spanbound=',best.fit$error.bound,', C=',best.fit$C,',', length(best.fit$features),'features.\n')
  
  result <- list(features=feat, error.bound = best.bound, fit=best.fit)
  class(result) <- 'rfe'
  return(result)
}

#' Predict Method for RFE Fits
#'
#' Obtains predictions from a fitted RFE object.
#'
#' @param object a fitted object of class inheriting from 'rfe'
#' @param newdata a matrix with variables to predict
#' @param type \code{response} gives the predictions \code{class} gives the predicted classes.
#' @param ... currently ignored.
#' @return the predictions.
#' @export
#' @callGraphPrimitives
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' \dontrun{
#' library(pathClass)
#' data(example_data)
#' fit = fit.rfe(x[1:5,], y[1:5], DEBUG=T)
#' predict(fit, newdata=x[6:10,])
#' }
predict.rfe = function(object, newdata, type="response", ...){
    svm.predict(object$fit, newdata, type)
}


