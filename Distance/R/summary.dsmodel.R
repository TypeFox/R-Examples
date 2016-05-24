#' Summary of distance sampling analysis 
#' 
#' Provides a brief summary of a distance sampling analysis. This includes 
# parameters, model selection criterion, and optionally abundance in the
# covered (sampled) region and its standard error.
# 
# The argument \code{N} is used to suppress computation of
# abundance and average detection probability in calls to summarize the
# \code{ds} and either \code{io.fi} or \code{trial.fi} for summaries of
# \code{io} and \code{trial} objects respectively which are composed of a
# \code{ds} model object and a mark-recapture model object. The corresponding
# print function is called to print the summary results.
# 
#' @aliases summary.dsmodel
#' @param object a distance analysis
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return list of extracted and summarized objects
#' @note This function just calls \code{\link{summary.ds}} and 
#'       \code{\link{dht}}, collates and prints the results in a nice
#'       way.
#' @author David L. Miller
#' @keywords utility
#' @export
summary.dsmodel <- function(object,...){

  #  se if TRUE, computes standard errors
  #  N if TRUE, computes abundance in covered (sampled) region
  ans <- list(ds=summary(object$ddf,se=TRUE,N=TRUE), 
              dht=object$dht, 
              ddf=object$ddf)

  class(ans) <- "summary.dsmodel"
  return(ans)
}
