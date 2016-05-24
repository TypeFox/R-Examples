#' Plot the sensitivity, specificity, accuracy and roc curves.
#'
#' This function plots the (partial) sensitivity, specificity, accuracy and roc curves.
#' 
#' @param x an object produced by one of the functions \code{sensitivity}, \code{specificity},  \code{accuracy}, or \code{roc}
#' @param y Not used.
#' @param ... Arguments to be passed to methods, such as graphical parameters. See ?plot
#' @param type Type of plot. Default is line plot.
#' @param add Logical. If TRUE the curve is added to an existing plot. If FALSE a new plot is created.
#' @param min a numeric value between 0 and 1, denoting the cutoff that defines the start of the area under the curve
#' @param max a numeric value between 0 and 1, denoting the cutoff that defines the end of the area under the curve
#' 
#' @examples
#' 
#' data(churn)
#' 
#' plot(sensitivity(churn$predictions,churn$labels))
#' 
#' plot(specificity(churn$predictions,churn$labels))
#' 
#' plot(accuracy(churn$predictions,churn$labels))
#' 
#' plot(roc(churn$predictions,churn$labels))
#' 
#' 
#' @references Ballings, M., Van den Poel, D., Threshold Independent Performance Measures for Probabilistic Classifcation Algorithms, Forthcoming.
#' @seealso \code{\link{sensitivity}}, \code{\link{specificity}}, \code{\link{accuracy}}, \code{\link{roc}}, \code{\link{auc}}, \code{\link{plot}}
#' @author Authors: Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@UGent.be}
#' @method plot AUC
plot.AUC <- function(x,y=NULL, ...,type='l',add=FALSE, min=0,max=1) {
  
  if (any(class(x) == "roc")) {
    
    if (min != 0 || max != 1 ) {
      x$fpr <- x$fpr[x$cutoffs >= min & x$cutoffs <= max]
      x$tpr <- x$tpr[x$cutoffs >= min & x$cutoffs <= max]
    }
    
  }else{
    
    if (min != 0 || max != 1 ) {
      ind <- x$cutoffs >= min & x$cutoffs <= max
      x$cutoffs <- x$cutoffs[ind]
      x$measure <- x$measure[ind]
    }
  } 
  
  
  if (any(class(x) == "roc")) {
  
          if (add==FALSE) {
            plot(x$fpr,x$tpr, type=type,xlab='1- specificity', ylab='sensitivity',xlim=c(0,1), ylim=c(0,1),...)
            lines(x=seq(0,1,by=0.01),y=seq(0,1,by=0.01),lty=1, col='grey')
          }else {
            lines(x$fpr,x$tpr,type=type, xlab='1- specificity', ylab='sensitivity',...)  
          }
  }else if (any(class(x) == "accuracy")) {

          
          if (add==FALSE) {
            plot(x$cutoffs,x$measure,type=type, xlab='Cutoffs', ylab='Accuracy',xlim=c(0,1), ylim=c(0,1),...)
          }else {
            lines(x$cutoffs,x$measure,type=type, xlab='Cutoffs', ylab='Accuracy', ...)   
          }
  }else if (any(class(x) == "specificity")) {
  
          
          if (add==FALSE) {
            plot(x$cutoffs,x$measure,type=type, xlab='Cutoffs', ylab='Specificity',xlim=c(0,1), ylim=c(0,1),...)
          }else {
            lines(x$cutoffs,x$measure,type=type, xlab='Cutoffs', ylab='Specificity', ...)   
          }
  }else if (any(class(x) == "sensitivity")) {
          
          
          if (add==FALSE) {
            plot(x$cutoffs, x$measure,type=type, xlab='Cutoffs', ylab='Sensitivity',xlim=c(0,1), ylim=c(0,1), ...)
          }else {
            lines(x$cutoffs, x$measure,type=type, xlab='Cutoffs', ylab='Sensitivity', ...)   
          }
  }
}