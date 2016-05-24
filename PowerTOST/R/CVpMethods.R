#------------------------------------------------------------------------------
# print and summary S3 methods for nice printing of pooled CV
# 
# Author: dlabes
#------------------------------------------------------------------------------

print.CVp <- function(x, digits = 4, verbose=FALSE, ...)
{
  if (verbose){
    cat("Pooled CV = ",format(x$CV,digits=digits), 
        " with ", x$df, " degrees of freedom", sep="")
    if (!is.null(x$robust)) {if (x$robust) cat(" (robust df's)")}
    cat("\n")
    cat("Upper ",(1-x$alpha)*100,"% confidence limit of CV = ", 
        format(x$CVupper, digits=digits),"\n",sep="")
  } else {
    cat(format(x$CV,digits=digits), " with ", x$df, 
        " degrees of freedom", sep="")
    if (!is.null(x$robust)) {if (x$robust) cat(" (robust df's)")}
    cat("\n")
  }
  invisible(x)
}
