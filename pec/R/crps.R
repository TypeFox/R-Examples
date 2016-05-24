#' Summarizing prediction error curves
#' 
#' Computes the cumulative prediction error curves, aka integrated Brier
#' scores, in ranges of time.
#' 
#' The cumulative prediction error (continuous ranked probability score) is
#' defined as the area under the prediction error curve, hence the alias name,
#' ibs, which is short for integrated Brier score.
#' 
#' @aliases crps ibs
#' @param object An object with estimated prediction error curves obtained with
#' the function \link{pec}
#' @param models Which models in \code{object$models} should be considered.
#' @param what The name of the entry in \code{x} to be cumulated. Defauls to
#' \code{PredErr} Other choices are \code{AppErr}, \code{BootCvErr},
#' \code{Boot632}, \code{Boot632plus}.
#' @param times Time points at which the integration of the prediction error
#' curve stops.
#' @param start The time point at which the integration of the prediction error
#' curve is started.
#' @return A matrix with a column for the crps (ibs) at every requested time
#' point and a row for each model
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{pec}}
#' @references E. Graf et al.  (1999), Assessment and comparison of prognostic
#' classification schemes for survival data. Statistics in Medicine, vol 18,
#' pp= 2529--2545.
#' 
#' Gerds TA, Cai T & Schumacher M (2008) The performance of risk prediction
#' models Biometrical Journal, 50(4), 457--479
#' @keywords survival
#' @examples
#' set.seed(18713)
#' library(prodlim)
#' library(survival)
#' dat=SimSurv(100)
#' pmodel=coxph(Surv(time,status)~X1+X2,data=dat)
#' perror=pec(list(Cox=pmodel),Hist(time,status)~1,data=dat)
#' 
#' ## cumulative prediction error
#' crps(perror,times=1) # between min time and 1
#' ## same thing:
#' ibs(perror,times=1) # between min time and 1
#' crps(perror,times=1,start=0) # between 0 and 1
#' crps(perror,times=seq(0,1,.2),start=0) # between 0 and seq(0,1,.2)
#' 
#' 
#' @export
crps <- function(object,
                 models,
                 what,
                 times,
                 start){
  stopifnot(class(object)[1] == "pec")
  
  # {{{find the prediction models
  if (missing(models)) models <- 1:length(object$models)
  else
    if (!is.numeric(models))
      models <- names(object$models)[match(models,names(object$models))]
  # }}}
  # {{{times
  object.times <- object$time
  if(missing(times)) times <- object$maxtime
  if (any(times>object$maxtime)) {
    warning(paste("You asked to integrate until times where prediction error curves are not defined.", format(object$maxtime,nsmall=2,digits=2)))
    times <- times[times<=object$maxtime]
  }
  ## if (!(object$exact))
  ## warning("Exact Only ", length(object.times)," time point",ifelse(length(times)==1,"","s")," used for computation of ")
  ##  time range
  if (missing(start)) start <- object$start
  # }}}
  # {{{ what errors
  if (missing(what) || is.null(what)){
    what <- grep(c("Err$"),names(object),value=TRUE)
  }
  # }}}
  # {{{ for each element of what: evaluate crps at times
  out <- lapply(what,function(w){
    est <- object[[w]][models]
    y <- sapply(times,function(t){
      intx <- sapply(est, function(y){
        Dint(x=object.times,
             y=y,
             range=c(start,t))
      })
    })
    if (!is.null(dim(y))){
      tnames <- paste("time=",round(times,1),sep="")
      tnames[times<1] <- paste("time=",signif(times[times<1],2),sep="")
      colnames(y) <- paste("IBS[",start,";",tnames,")",sep="")
      y}
    else{
      y
    }
  })
  # }}}
  # {{{ prepare output
  NW <- length(what)
  NT <- length(times)
  if (NW==1)
    out <- out[[1]]
  else
    names(out) <- what
  if (NT==1){
    if(NW>1){
      out <- do.call("cbind",out)
      colnames(out) <- what
    }
  }
  # }}}
  class(out) <- "crps"
  out
}

