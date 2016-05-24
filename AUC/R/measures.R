#' Compute the receiver operating characteristic (ROC) curve.
#'
#' This function computes the receiver operating characteristic (ROC) curve required for the \code{auc} function and the \code{plot} function.
#' 
#' @param predictions A numeric vector of classification probabilities (confidences, scores) of the positive event.
#' @param labels A factor of observed class labels (responses) with the only allowed values \{0,1\}.
#' 
#' @examples
#' 
#' data(churn)
#' 
#' roc(churn$predictions,churn$labels)
#' 
#' @references Ballings, M., Van den Poel, D., Threshold Independent Performance Measures for Probabilistic Classifcation Algorithms, Forthcoming.
#' @seealso \code{\link{sensitivity}}, \code{\link{specificity}}, \code{\link{accuracy}}, \code{\link{roc}}, \code{\link{auc}}, \code{\link{plot}}
#' @return A list containing the following elements:
#' \item{cutoffs}{A numeric vector of threshold values}
#' \item{fpr}{A numeric vector of false positive rates corresponding to the threshold values}
#' \item{tpr}{A numeric vector of true positive rates corresponding to the threshold values}
#' @author Authors: Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@UGent.be}
roc <-  function(predictions, labels) {
                
                #This function is a scaled down faster version of parts of code of the ROCR package.
                #It has less functionality and less error handling and is focused on speed.
                #please see ROCR
  
                cm <- .confusionMatrix(predictions,labels,FALSE,'roc')
                
                x <- cm$fp / cm$n.neg
                y <- cm$tp / cm$n.pos
                
                finite.bool <- is.finite(x) & is.finite(y)
                x <- x[ finite.bool ]
                y <- y[ finite.bool ]
                if (length(x) < 2) {
                  stop(paste("Not enough distinct predictions to compute area",
                             "under the ROC curve."))
                } 
                ans <- list(cutoffs=cm$cutoffs, fpr=x, tpr=y )
                class(ans) <- c('AUC','roc')
                return(ans)  
                
}
            



#' Compute the sensitivity curve.
#'
#' This function computes the sensitivity curve required for the \code{auc} function and the \code{plot} function.
#' 
#' @param predictions A numeric vector of classification probabilities (confidences, scores) of the positive event.
#' @param labels A factor of observed class labels (responses) with the only allowed values \{0,1\}.
#' @param perc.rank A logical. If TRUE (default) the percentile rank of the predictions is used.
#' 
#' @examples
#' 
#' data(churn)
#' 
#' sensitivity(churn$predictions,churn$labels)
#' 
#' @references Ballings, M., Van den Poel, D., Threshold Independent Performance Measures for Probabilistic Classifcation Algorithms, Forthcoming.
#' @seealso \code{\link{sensitivity}}, \code{\link{specificity}}, \code{\link{accuracy}}, \code{\link{roc}}, \code{\link{auc}}, \code{\link{plot}}
#' @return A list containing the following elements:
#' \item{cutoffs}{A numeric vector of threshold values}
#' \item{measure}{A numeric vector of sensitivity values corresponding to the threshold values}
#' @author Authors: Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@UGent.be}
sensitivity <-  function(predictions, labels, perc.rank=TRUE) {
                
                cm <- .confusionMatrix(predictions,labels,perc.rank, 'sensitivity')
                
                ans <- list( cutoffs=c(cm$cutoffs,0), measure=c(cm$tp / cm$n.pos,1) )
                class(ans) <- c('AUC','sensitivity')
                return(ans)
}
            



 
#' Compute the specificity curve.
#'
#' This function computes the specificity curve required for the \code{auc} function and the \code{plot} function.
#' 
#' @param predictions A numeric vector of classification probabilities (confidences, scores) of the positive event.
#' @param labels A factor of observed class labels (responses) with the only allowed values \{0,1\}.
#' @param perc.rank A logical. If TRUE (default) the percentile rank of the predictions is used.
#' 
#' @examples
#' 
#' data(churn)
#' 
#' specificity(churn$predictions,churn$labels)
#' 
#' @references Ballings, M., Van den Poel, D., Threshold Independent Performance Measures for Probabilistic Classifcation Algorithms, Forthcoming.
#' @seealso \code{\link{sensitivity}}, \code{\link{specificity}}, \code{\link{accuracy}}, \code{\link{roc}}, \code{\link{auc}}, \code{\link{plot}}
#' @return A list containing the following elements:
#' \item{cutoffs}{A numeric vector of threshold values}
#' \item{measure}{A numeric vector of specificity values corresponding to the threshold values}
#' @author Authors: Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@UGent.be}
specificity <-  function(predictions, labels, perc.rank=TRUE) {
                
                cm <- .confusionMatrix(predictions,labels, perc.rank, 'specificity')
                
                ans <- list( cutoffs=c(cm$cutoffs,0), measure=c(cm$tn / cm$n.neg,0) )
                
                class(ans) <- c('AUC','specificity')
                return(ans)
}




#' Compute the accuracy curve.
#'
#' This function computes the accuracy curve required for the \code{auc} function and the \code{plot} function.
#' 
#' @param predictions A numeric vector of classification probabilities (confidences, scores) of the positive event.
#' @param labels A factor of observed class labels (responses) with the only allowed values \{0,1\}.
#' @param perc.rank A logical. If TRUE (default) the percentile rank of the predictions is used.
#' 
#' @examples
#' 
#' data(churn)
#' 
#' accuracy(churn$predictions,churn$labels)
#' 
#' @references Ballings, M., Van den Poel, D., Threshold Independent Performance Measures for Probabilistic Classifcation Algorithms, Forthcoming.
#' @seealso \code{\link{sensitivity}}, \code{\link{specificity}}, \code{\link{accuracy}}, \code{\link{roc}}, \code{\link{auc}}, \code{\link{plot}}
#' @return A list containing the following elements:
#' \item{cutoffs}{A numeric vector of threshold values}
#' \item{measure}{A numeric vector of accuracy values corresponding to the threshold values}
#' @author Authors: Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@UGent.be}
accuracy <-  function(predictions, labels, perc.rank=TRUE) {
                
                cm <- .confusionMatrix(predictions,labels, perc.rank, 'accuracy')
                
                #at cutoff of 0 the accuracy equals maximal tpr (i.e., 1) times the proportion of positives
                ans <- list( cutoffs=c(cm$cutoffs,0), measure= c((cm$tn+cm$tp) / (cm$n.pos + cm$n.neg), mean(as.integer(as.character(labels)) ) ))
                
                class(ans) <- c('AUC','accuracy')
                return(ans)
}