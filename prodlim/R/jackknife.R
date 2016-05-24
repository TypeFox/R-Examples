#' Compute jackknife pseudo values.
#' 
#' Compute jackknife pseudo values based on marginal Kaplan-Meier estimate of
#' survival, or based on marginal Aalen-Johansen estimate of cumulative
#' incidence.
#' 
#' @title Compute jackknife pseudo values.
#' @aliases jackknife jackknife.survival jackknife.competing.risks
#' @param object Object of class \code{"prodlim"}.
#' @param times Time points at which to compute pseudo values.
#' @param cause For competing risks the cause of failure.
#' @param keepResponse If \code{TRUE} add the model response,
#' i.e. event time, event status, etc. to the result.
#' @param ... not used
#' @note The R-package pseudo does a similar job, and appears to be a little faster in small samples, but much slower in large samples. See examples.
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{prodlim}}
#' @references Andersen PK & Perme MP (2010). Pseudo-observations in survival
#' analysis Statistical Methods in Medical Research, 19(1), 71-99.
#' @keywords survival
##' @examples
##' 
##' 
##' ## pseudo-values for survival models
##' 
##' d=SimSurv(20) 
##' f=prodlim(Hist(time,status)~1,data=d) 
##' jackknife(f,times=c(3,5))
##' 
##' ## in some situations it may be useful to attach the
##' ## the event time history to the result
##' jackknife(f,times=c(3,5),keepResponse=TRUE)
##' 
##' # pseudo-values for competing risk models
##' d=SimCompRisk(10) 
##' f=prodlim(Hist(time,event)~1,data=d) 
##' jackknife(f,times=c(3,10),cause=1)
##' jackknife(f,times=c(3,10,17),cause=2)
##' 
#' @export
jackknife <- function(object,times,cause,keepResponse=FALSE,...){
  if (object$model=="survival")
    jackknife.survival(object=object,times=times,keepResponse=keepResponse,...)
  else if (object$model=="competing.risks")
    jackknife.competing.risks(object=object,
                              times=times,
                              cause=cause,
                              keepResponse=keepResponse,
                              ...)
  else stop("No method for jackknifing this object.")
}

#' @export
jackknife.survival <- function(object,times,keepResponse=FALSE,...){
  S <- predict(object,times=times,newdata=object$model.response)
  Sk <- leaveOneOut.survival(object,times,...)
  N <- NROW(Sk)
  Jk <- t(N*S-t((N-1)*Sk))
  colnames(Jk) <- paste("t",times,sep=".")
  if (keepResponse==TRUE){
    Jk <- cbind(object$model.response,Jk)
  }
  ## re-order the pseudo-values  
  Jk <- Jk[object$originalDataOrder,,drop=FALSE]
  Jk
}
#' @export
jackknife.competing.risks <- function(object,times,cause,keepResponse=FALSE,...){
  F <- predict(object,times=times,newdata=object$model.response,cause=cause)
  Fk <- leaveOneOut.competing.risks(object,times,cause=cause,...)
  N <- NROW(Fk)
  Jk <- t(N*F-t((N-1)*Fk))
  colnames(Jk) <- paste("t",times,sep=".")
  if (keepResponse==TRUE){
    Jk <- cbind(object$model.response,Jk)
    colnames(Jk)[(NCOL(Jk)-length(times)+1):NCOL(Jk)] <- paste("t",times,sep=".")
  }
  ## re-order the pseudo-values
  Jk <- Jk[object$originalDataOrder,,drop=FALSE]
  Jk
}



