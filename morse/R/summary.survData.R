#' Summary for survData objects
#' 
#' The generic \code{summary} S3 method for the \code{survData} class provides
#' information about the structure of the dataset and the experimental design.
#' 
#' @param object an object of class \code{survData}
#' @param quiet when \code{FALSE}, prints summary on standard output
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @return The function returns a list with the following fields:
#' \item{NbrepTimeConc}{nb of replicates for all concentrations and time points}
#' \item{NbsurvTimeConc}{nb of surviving ind. for all concentrations and time points}

#' @seealso \code{\link{survData}}
#' 
#' @examples
#' # (1) Load the data
#' data(cadmium1)
#' 
#' # (2) Create a survData object
#' cadmium1 <- survData(cadmium1)
#' 
#' # (3) Summarize the dataset
#' summary(cadmium1)
#' 
#' @keywords summary
#' 
#' @export
summary.survData <- function(object, quiet = FALSE, ...) {
  # matrix of number of replicate by time / conc
  ans1 <- table(object[, c("conc", "time")])
  
  # matrix of number of survival (sum of all replicate) by time / conc
  ans2 <- tapply(object$Nsurv, list(as.factor(object$conc),
                                    as.factor(object$time)), sum)
  
  if (! quiet) {
    cat("\nNumber of replicate per time and concentration: \n")
    print(ans1)
    cat("\nNumber of survival (sum of replicate) per time and concentration: \n")
    print(ans2)
  }
  
  invisible(list(NbrepTimeConc = ans1,
                 NbsurvTimeConc = ans2))
}
