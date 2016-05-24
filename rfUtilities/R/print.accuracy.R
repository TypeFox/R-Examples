#' @title Print accuracy
#' @description print method for class "accuracy"
#' @param x    Object of class accuracy
#' @param ...  Ignored
#'
#' @method print accuracy
#'
#' @export
print.accuracy <- function(x, ...) {
  if( dim(x$confusion)[1] > 2) {   
    cat("Accuracy (PCC):", paste(x$PCC, "%", sep=""), "\n\n")
    cat("Cohen's Kappa:", x$kappa, "\n\n")
    cat("Users accuracy:", "\n")
      print(x$users.accuracy)
        cat("", "\n\n")
    cat("Producers accuracy:", "\n")
      print(x$producers.accuracy)
        cat("", "\n\n")
    cat("Confusion matrix", "\n")
      print(x$confusion)
    } else {
    cat("Accuracy (PCC):", paste(x$PCC, "%", sep=""), "\n\n")
    cat("Cohen's Kappa:", x$kappa, "\n\n")
	cat("Area under the ROC curve:", x$auc, "\n\n")
    cat("Users accuracy:", "\n")
      print(x$users.accuracy)
        cat("", "\n\n")
    cat("Producers accuracy:", "\n")
      print(x$producers.accuracy)
        cat("", "\n\n")
	cat("True Skill statistic:", x$true.skill, "\n\n")
    cat("Sensitivity:", x$sensitivity, "\n\n")
    cat("Sensitivity:", x$specificity, "\n\n")
    cat("Positive Likelihood Ratio:", x$plr, "\n\n")
    cat("Negative Likelihood Ratio:", x$nlr, "\n\n")
    cat("Type I error:", x$typeI.error, "\n\n")
    cat("Type II error:", x$typeII.error, "\n\n")
    cat("F-score:", x$f.score, "\n\n")
	cat("Matthews correlation coefficient:", x$matthews, "\n\n")
    cat("Confusion matrix", "\n")
      print(x$confusion)	
  }	
}  
