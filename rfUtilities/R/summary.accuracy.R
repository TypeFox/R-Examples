#' @title Summarizing accuracy
#' @description Summary method for class "accuracy".
#' @param object  Object of class accuracy
#' @param ... Ignored
#'
#' @method summary accuracy
#'
#' @export
summary.accuracy <- function(object, ...) {
  if( dim(object$confusion)[1] > 2) {   
    cat("Accuracy (PCC):", paste(object$PCC, "%", sep=""), "\n\n")
    cat("Cohen's Kappa:", object$kappa, "\n\n")
    cat("Users accuracy:", "\n")
      print(object$users.accuracy)
        cat("", "\n\n")
    cat("Producers accuracy:", "\n")
      print(object$producers.accuracy)
        cat("", "\n\n")
    cat("Confusion matrix", "\n")
      print(object$confusion)
    } else {
    cat("Accuracy (PCC):", paste(object$PCC, "%", sep=""), "\n\n")
    cat("Cohen's Kappa:", object$kappa, "\n\n")
	cat("Area under the ROC curve:", object$auc, "\n\n")
    cat("Users accuracy:", "\n")
      print(object$users.accuracy)
        cat("", "\n\n")
    cat("Producers accuracy:", "\n")
      print(object$producers.accuracy)
        cat("", "\n\n")
    cat("True Skill statistic:", object$true.skill, "\n\n")
    cat("Sensitivity:", object$sensitivity, "\n\n")
    cat("Sensitivity:", object$specificity, "\n\n")
    cat("Positive Likelihood Ratio:", object$plr, "\n\n")
    cat("Negative Likelihood Ratio:", object$nlr, "\n\n")
    cat("Type I error:", object$typeI.error, "\n\n")
    cat("Type II error:", object$typeII.error, "\n\n")
    cat("F-score:", object$f.score, "\n\n")
	cat("Matthews correlation coefficient:", object$matthews, "\n\n")
    cat("Confusion matrix", "\n")
      print(object$confusion)	
  }	
}
