#' Print summary of distance detection function model object
#'
#' Provides a brief summary of data and fitted detection probability model
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error. What is printed depends
#' on the corresponding call to summary.
#'
#' @aliases print.summary.io
#' @export
#' @param x a summary of \code{ddf} model object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @author Jeff Laake
#' @seealso \code{\link{summary.io}}
#' @keywords utility
print.summary.io <- function(x,...){

  print(x$mr.summary)
  cat("\n")
  print(x$ds.summary)
  cat("\n\nSummary for io object\n")
  cat("Total AIC value : ",x$AIC,"\n\n")

  parameters <- data.frame(Estimate=c(x$average.p,x$Nhat))
  row.names(parameters) <- c("Average p", "N in covered region")
  if(!is.null(x$average.p.se)){
    parameters$SE <- c(x$average.p.se,x$Nhat.se)
    parameters$CV <- parameters$SE/parameters$Estimate
  }
  print(parameters)
  invisible(NULL)
}
