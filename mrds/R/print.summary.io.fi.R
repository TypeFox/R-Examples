#' Print summary of distance detection function model object
#'
#' Provides a brief summary of data and fitted detection probability model
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error. What is printed depends
#' on the corresponding call to summary.
#'
#' @export
#' @param x a summary of \code{ddf} model object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @author Jeff Laake
#' @seealso \code{\link{summary.io.fi}}
#' @keywords utility
print.summary.io.fi <- function(x,...){

  cat("\nSummary for io.fi object \n")
  cat("Number of observations   : ", x$n,"\n")
  cat("Number seen by primary   : ", x$n1,"\n")
  cat("Number seen by secondary : ", x$n2,"\n")
  cat("Number seen by both      : ", x$n3,"\n")
  cat("AIC                      : ", x$aic, "\n")

  cat("\n\nConditional detection function parameters:\n")
  print(x$cond.det.coef)

  cat("\n")
  if(!is.null(x$Nhat)){
    parameters <- data.frame(Estimate=c(x$average.p,x$average.p0.1,
                                        x$average.p0.2,x$average.p0,x$Nhat))
    row.names(parameters) <- c("Average p","Average primary p(0)",
                               "Average secondary p(0)","Average combined p(0)",
                               "N in covered region")
    if(!is.null(x$average.p.se)){
      parameters$SE <- c(x$average.p.se,x$average.p0.1.se,x$average.p0.2.se,
                           x$average.p0.se,x$Nhat.se)
      parameters$CV <- parameters$SE/parameters$Estimate
    }
  }else{
    parameters <- data.frame(Estimate=c(x$average.p0.1,x$average.p0.2,
                           x$average.p0))
    row.names(parameters) <- c("Average primary p(0)", "Average secondary p(0)",
                            "Average combined p(0)")
    if(!is.null(x$average.p0.se)){
      parameters$SE <- c(x$average.p0.1.se,x$average.p0.2.se,
                           x$average.p0.se)
      parameters$CV <- parameters$SE/parameters$Estimate
    }
  }
  print(parameters)
  invisible(NULL)
}
