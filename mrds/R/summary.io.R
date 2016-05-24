#' Summary of distance detection function model object
#'
#' Provides a brief summary of data and fitted detection probability model
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error.
#'
#' The argument \code{N} is used to suppress computation of
#' abundance and average detection probability in calls to summarize the
#' \code{ds} and either \code{io.fi} or \code{trial.fi} for summaries of
#' \code{io} and \code{trial} objects respectively which are composed of a
#' \code{ds} model object and a mark-recapture model object. The corresponding
#' print function is called to print the summary results.
#'
#' @export
#' @param object a \code{ddf} model object
#' @param se if TRUE, computes standard errors
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return list of extracted and summarized objects
#' @note This function is called by the generic function \code{summary} for any
#'   \code{ddf} model object.  Each function can be called directly by the
#'   user, but it is typically safest to use the generic function
#'   \code{summary} which calls the appropriate function based on the type of
#'   \code{ddf} model.
#' @author Jeff Laake
#' @keywords utility
summary.io <- function(object,se=TRUE,...){
  # Uses: predict.io, summary.ds, summary.io.fi

  model <- object
  avgp <- function(model,pdot,...){return(pdot)}
  ddfobj <- model$ds$ds$aux$ddfobj
  n <- nrow(ddfobj$xmat)
  ans <- list(mr.summary  = summary(model$mr,se=se,N=FALSE,model,ddfobj),#io.fi
              ds.summary  = summary(model$ds,se=se,N=FALSE), #ds
              Nhat        = model$Nhat,
              AIC         = model$criterion,
              average.p   = n/model$Nhat,
              mono        = model$ds$aux$mono,
              mono.strict = model$ds$aux$mono.strict)

  # calculate standard errors for p and Nhat (covered area)
  if(se){
    se.obj <- calc.se.Np(model, avgp, n, ans$average.p)

    ans$average.p.se <- se.obj$average.p.se
    ans$Nhat.se <- se.obj$Nhat.se

    ans$cv <- ans$Nhat.se/model$Nhat
  }
  class(ans) <- "summary.io"
  return(ans)
}
