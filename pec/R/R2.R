#' Explained variation for survival models
#' 
#' This function computes a time-dependent $R^2$ like measure of the variation
#' explained by a survival prediction model, by dividing the mean squared error
#' (Brier score) of the model by the mean squared error (Brier score) of a
#' reference model which ignores all the covariates.
#' 
#' 
#' In survival analysis the prediction error of the Kaplan-Meier estimator
#' plays a similar role as the total sum of squares in linear regression.
#' Hence, it is a sensible reference model for $R^2$.
#' 
#' @param object An object with estimated prediction error curves obtained with
#' the function \link{pec}
#' @param models For which of the models in \code{object$models} should we
#' compute $R^2(t). By default all models are used except for the reference
#' model.
#' @param what The name of the entry in \code{x} to be used. Defauls to
#' \code{PredErr} Other choices are \code{AppErr}, \code{BootCvErr},
#' \code{Boot632}, \code{Boot632plus}.
#' @param times Time points at which the summaries are shown.
#' @param reference Position of the model whose prediction error is used as the
#' reference in the denominator when constructing $R^2$
#' @return A matrix where the first column holds the times and the following
#' columns are the corresponding $R^2$ values for the requested prediction
#' models.
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
#' 
#' set.seed(18713)
#' library(prodlim)
#' library(survival)
#' dat=SimSurv(100)
#' nullmodel=prodlim(Hist(time,status)~1,data=dat)
#' pmodel1=coxph(Surv(time,status)~X1+X2,data=dat)
#' pmodel2=coxph(Surv(time,status)~X2,data=dat)
#' perror=pec(list(Cox1=pmodel1,Cox2=pmodel2),Hist(time,status)~1,data=dat,reference=TRUE)
#' R2(perror,times=seq(0,1,.1),reference=1)
#' 
#' @export
R2 <- function(object,
               models,
               what,
               times,
               reference=1){
  
  stopifnot(class(object)[1] == "pec")
  
  # {{{find the prediction models
  
  if (missing(models))
    models <- (1:length(object$models))[-reference]
  else
    if (!is.numeric(models))
      models <- match(models,names(object$models))
  # }}}
  # {{{ what errors
  if (missing(what) || is.null(what)){
    what <- grep(c("Err$"),names(object),value=TRUE)
  }
  # }}}
  # {{{ find the times

  object.times <- object$time
  if(missing(times)) times <- object$maxtime
  if (!(object$exact || length(object.times)>100))
    warning("Only ", length(time)," time point",ifelse(length(times)==1,"","s")," used")
  # }}}
  # {{{ for each element of what: evaluate R2 at times

  out <- lapply(what,function(e){
    if (is.null(object[[e]])) stop("No values for computing R^2")
    ref.error <- object[[e]][[reference]]
    out <- data.frame(do.call("cbind",lapply(1:length(models),function(w){
      rr <- 1-object[[e]][[models[w]]]/ref.error
      rr[ref.error==0] <- 0
      rr
    })))
    names(out) <- names(object$models)[models]
    ## cat("R^2 based on the estimate stored in ",what,":\n\n")
    ## print(cbind(time=times,RR=rbind(0,out)[1+prodlim::sindex(object.times,times),,drop=FALSE]))
    cbind(time=times,RR=rbind(0,out)[1+prodlim::sindex(object.times,times),,drop=FALSE])
  })
  # }}}
  # {{{ prepare output
  NW <- length(what)
  NT <- length(times)
  names(out) <- what
  # }}}
  attr(out,"reference") <- names(object$models)[reference]
  class(out) <- "R2"
  out
}
