#' Print of \code{reproFitTT} object
#' 
#' This is the generic \code{print} S3 method for the \code{survFitTT} class.
#' It prints the underlying JAGS model and some information on the Bayesian 
#' inference procedure.
#' 
#' @param x An object of class \code{reproFitTT}
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @seealso \code{\link{reproFitTT}}
#' 
#' @examples
#' # (1) Load the data
#' data(cadmium1)
#' 
#' # (2) Create a reproData object
#' cadmium1 <- reproData(cadmium1)
#' 
#' \dontrun{
#' # (3) Run the reproFitTT function with the log-logistic
#' # model
#' out <- reproFitTT(cadmium1, ecx = c(5, 10, 15, 20, 30, 50, 80),
#' quiet = TRUE)
#' 
#' # (4) Print the reproFitTT object
#' out
#' }
#' 
#' @keywords print
#' 
#' @export
print.reproFitTT <- function(x, ...) {
  # print the model text and the Jags Computing information
  # for an object of class reproFitTT
  
  # M.C.M.C. informations
  cat("Model:\n")
  print(x$model)
  cat("\nComputing information:\n\n")
  cat("\n", "Iterations = ", x$n.iter[["start"]], ":",
      x$n.iter[["end"]], "\n", sep = "")
  cat("Thinning interval =", x$n.thin, "\n")
  cat("Number of chains =", x$n.chains, "\n")
  cat("Sample size per chain =",
      (x$n.iter[["end"]] - x$n.iter[["start"]]) / x$n.thin + 1, "\n")
}
