#' @title Print random forests cross-validation
#' @description Print method for rf.cv objects
#'    
#' @param x    Object of class rf.cv
#' @param ...  Ignored
#'
#' @method print rf.cv
#'
#' @export
print.rf.cv <- function(x, ...) {
  cat("Classification accuracy for cross-validation", "\n")
  cv <- data.frame()	
    cv <- rbind(cv,
            apply(x$cross.validation$cv.users.accuracy, MARGIN = 2, stats::median),
            apply(x$cross.validation$cv.producers.accuracy, MARGIN = 2, stats::median))
		      row.names(cv) <- c("users.accuracy", "producers.accuracy")
              names(cv) <- names(x$cross.validation$cv.users.accuracy)
    cat("", "\n")
	print( cv )
	cat("", "\n")
	cat("Cross-validation Kappa", "=", stats::median(x$cross.validation$cv.oob[,"kappa"]), "\n")
	cat("Cross-validation OOB Error", "=", stats::median(x$cross.validation$cv.oob[,"OOB"]), "\n")
	cat("Cross-validation error variance", "=", stats::var(x$cross.validation$cv.oob[,"OOB"]), "\n")
    cat("", "\n")
	cat("", "\n")
	
  cat("Classification accuracy for model", "\n")
  mdl <- data.frame()	
    mdl <- rbind(mdl,
            apply(x$model$model.users.accuracy, MARGIN = 2, stats::median),
            apply(x$model$model.producers.accuracy, MARGIN = 2, stats::median))
		      row.names(mdl) <- c("users.accuracy", "producers.accuracy")
              names(mdl) <- names(x$model$model.users.accuracy)
    cat("", "\n")
	print( mdl )
	cat("", "\n")
    cat("Model Kappa", "=", stats::median(x$model$model.oob[,"kappa"]), "\n")
	cat("Model OOB Error", "=", stats::median(x$model$model.oob[,"OOB"]), "\n")
	cat("Model error variance", "=", stats::var(x$model$model.oob[,"OOB"]), "\n")
}
