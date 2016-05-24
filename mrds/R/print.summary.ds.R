#' Print summary of distance detection function model object
#'
#' Provides a brief summary of data and fitted detection probability model
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error. What is printed depends
#' on the corresponding call to summary.
#'
#' @aliases print.summary.ds
#' @export
#' @param x a summary of \code{ddf} model object
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return NULL
#' @author Jeff Laake
#' @seealso \code{\link{summary.ds}}
#' @keywords utility
print.summary.ds <- function (x,...){
  cat("\nSummary for ds object \n")
  cat("Number of observations : ", x$n,"\n")
  cat("Distance range         : ", x$left, " - ",x$width,"\n")
  cat("AIC                    : ", x$aic, "\n")
  cat("\nDetection function:\n",model.description(x),"\n")
  cat("\nDetection function parameters", "\n")
  cat("Scale Coefficients: ", "\n")
  print(x$coeff$key.scale)

  if(x$key %in% c("gamma","hr","th1","th2")) {
    cat("\nShape parameters: ", "\n")
    print(x$coeff$key.shape)
  }

  if (!is.null(x$coeff$adj.parm)) {
     cat("\nAdjustment term parameter(s): ", "\n")
     print(x$coeff$adj.parm)
  }

  cat("\n")

  if(x$mono & x$mono.strict){
    cat("\nStrict monotonicity constraints were enforced.\n")
  }else if(x$mono){
    cat("\nMonotonicity constraints were enforced.\n")
  }

  if(!is.null(x$Nhat)){
    parameters <- data.frame(Estimate=c(x$average.p,x$Nhat))
    row.names(parameters) <- c("Average p", "N in covered region")
    if(!is.null(x$average.p.se)){
      parameters$SE <- c(x$average.p.se,x$Nhat.se)
      parameters$CV <- parameters$SE/parameters$Estimate
    }
  }else{
    parameters <- data.frame(Estimate=c(x$average.p))
    row.names(parameters) <- c("Average p")
    if(!is.null(x$average.p.se)){
      parameters$SE <- c(x$average.p.se)
      parameters$CV <- parameters$SE/parameters$Estimate
    }
  }

  print(parameters)

  # Remind the user that monotonicity constraints were enforced
  if(x$mono & x$mono.strict){
    cat("\nStrict monotonicity constraints were enforced.\n")
  }else if(x$mono){
    cat("\nMonotonicity constraints were enforced.\n")
  }

  invisible()
}

