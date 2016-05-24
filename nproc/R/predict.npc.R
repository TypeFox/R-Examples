#' Predicting the outcome of a set of new observations using the fitted npc
#' object.
#' @export
#' @param object fitted npc object using \code{npc}.
#' @param newx a set of new observations.
#' @param pred.score a vector of scores for the new observations. Used when method = 'custom'.
#' @param ... additional arguments.
#' @return A list containing the predicted label and score.
#' \item{pred.label}{Predicted label vector.}
#' \item{pred.score}{Predicted score vector.}
#' @seealso \code{\link{npc}} and \code{\link{nproc}}
#' @examples
#' n = 1000
#' x = matrix(rnorm(n*2),n,2)
#' c = 1+3*x[,1]
#' y = rbinom(n,1,1/(1+exp(-c)))
#' xtest = matrix(rnorm(n*2),n,2)
#' ctest = 1+3*xtest[,1]
#' ytest = rbinom(n,1,1/(1+exp(-ctest)))
#'
#' ##Use logistic classifier and the default type I error control with alpha=0.05
#' #fit = npc(x, y, method = 'logistic')
#' #pred = predict(fit,xtest)
#' #fit.score = predict(fit,x)
#' #accuracy = mean(pred$pred.label==ytest)
#' #cat('Overall Accuracy: ',  accuracy,'\n')
#' #ind0 = which(ytest==0)
#' #ind1 = which(ytest==1)
#' #typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
#' #cat('Type I error: ', typeI, '\n')
#' #typeII = mean(pred$pred.label[ind1]!=ytest[ind1]) #type II error on test set
#' #cat('Type II error: ', typeII, '\n')

predict.npc <- function(object, newx = NULL, pred.score = NULL, ...) {

  if (object$method == "custom") {
    if (is.null(pred.score)) {
      stop("pred.score needed for method \"custom\".")
    }
    pred.score = pred.score
    if (object$sign == TRUE) {
      pred.label = outer(pred.score, object$cutoff, ">")
    } else {
      pred.label = outer(pred.score, object$cutoff, "<=")
    }
  } else {
    colnames(newx) <- paste("x", 1:ncol(newx), sep = "")
    if(object$split<=1){
      pred = pred.npc.core(object,newx)
      pred.score = pred$pred.score
      pred.label = pred$pred.label
    } else {
      pred.label = pred.npc.core(object[[1]],newx)$pred.label

      for(i in 2:object$split){
        pred.label = pred.label + pred.npc.core(object[[i]],newx)$pred.label
      }
      pred.label = (pred.label/object$split>0.5)
    }
  }


  return(list(pred.label = pred.label, pred.score = pred.score))

}
