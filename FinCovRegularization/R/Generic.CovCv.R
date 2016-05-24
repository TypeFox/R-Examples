#' @title Display a useful description of a CovCv object
#'
#' @description
#' Display a useful description of a CovCv object
#'
#' @param object CovCv object to summarise
#' @param ... other arguments ignored (for compatibility with generic)
#' @keywords internal
#' @export

summary.CovCv <- function(object, ...) {
  cat(object$n.cv, "times Cross-Validation repeated for",
      object$regularization, "parameter selection \n")
  cat("Optimal Parameter Value:", object$parameter.opt, "\n")
  if (object$regularization %in% c("Hard Thresholding", "Soft Thresholding")) {
    cat("Thresholding Method: ",
        paste(object$method, "thresholding", sep=" "), "\n")
  }
  if (object$norm %in% c("F","f")) {
    cat("Measurement: Frobenius norm \n")
  } else if (object$norm %in% c("O","o")) {
    cat("Measurement: operator norm \n")
  }
}


#' @title plot CovCv object
#'
#' @description
#' Visualizes the results of covariance matrix regularization parameter tuning
#'
#' @param x CovCv object to plot
#' @param ... other arguments ignored (for compatibility with generic)
#' @keywords internal
#' @importFrom graphics plot
#' @export

plot.CovCv <- function(x, ...) {
  if (x$regularization == "Banding") {
    plot(x = 0:(length(x$cv.error)-1), y = x$cv.error,
         type = "l", xlab = "k.grid", ylab = "cv.error")
  } else if (x$regularization == "Tapering") {
    plot(x = 0:(length(x$cv.error)-1), y = x$cv.error,
         type = "l", xlab = "l.grid", ylab = "cv.error")
  } else if (x$regularization %in% c("Hard Thresholding", "Soft Thresholding")) {
    plot(x = x$threshold.grid, y = x$cv.error,
         type = "l", xlab = "threshold.grid", ylab = "cv.error")
  }
}


#' @title print CovCv object
#'
#' @description
#' print selected optimal parameter
#'
#' @param x CovCv object to plot
#' @param ... other arguments ignored (for compatibility with generic)
#' @keywords internal
#' @export

print.CovCv <- function(x, ...) {
  print(x$parameter.opt)
}



