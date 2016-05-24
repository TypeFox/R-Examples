#' Summary for reproData objects
#' 
#' The generic \code{summary} S3 method for the \code{reproData} class provides
#' information about the structure of the dataset and the experimental design.
#' 
#' @param object an object of class \code{reproData}
#' @param quiet if \code{TRUE}, does no prints
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @return The function returns a list with the same fields than 
#' \code{\link{summary.survData}} plus an additional one:
#' \item{NboffTimeConc}{nb of offspring for all concentrations and time points}
#' 
#' @seealso \code{\link{reproData}}, \code{\link{summary.survData}}
#' 
#' @examples
#' # (1) Load the data
#' data(cadmium1)
#' 
#' # (2) Create a reproData object
#' cadmium1 <- reproData(cadmium1)
#' 
#' # (3) Summarize the dataset
#' summary(cadmium1)
#' 
#' @keywords summary
#' 
#' @export
summary.reproData <- function(object, quiet = FALSE, ...) {
  res <- summary.survData(object, quiet = quiet)
  
  # matrix of number of offspring (sum of all replicate) by time / conc
  ans3 <- tapply(object$Nrepro,
                 list(as.factor(object$conc), as.factor(object$time)), sum)
  
  if (! quiet) {
    cat("\nNumber of offspring (sum of replicate) per time and concentration: \n")
    print(ans3)
  }
  
  invisible(c(res, list(NboffTimeConc = ans3)))
}
