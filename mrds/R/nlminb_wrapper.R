#' Wrapper around \code{nlminb}
#'
#' This is a wrapper around nlminb to use scaling, as this is not available (nor will it be) in \code{\link{optimx}}.
#'
#' @import optimx
#' @inheritParams optimx::nlminb
#'
#' @param par starting parameters
#' @param ll log likelihood function
#' @param ugr gradent function
#' @param lower lower bounds on parameters
#' @param upper upper bounds on parameters
#' @param mcontrol control options
#' @param hess hessian function
#' @param ddfobj detection function specification object
#' @param data the data
#' @param \dots anything else to pass to \code{ll}
#' @return \code{optimx} object
#'
#' @importFrom stats nlminb
#' @author David L Miller, modified from \code{optimx.run} by JC Nash, R Varadhan, G Grothendieck.
nlminb_wrapper <- function(par, ll, ugr=NULL, lower=NULL, upper=NULL,
                           mcontrol, hess=NULL, ddfobj, data, ...){

  # need to do this from optimx.R
  optcfg <- optimx.setup(par, ll, ugr, hess, lower, upper,
                "nlminb", itnmax=NULL, hessian=FALSE, mcontrol, ...)
  mcontrol <- optcfg$ctrl

  ## most comments here are from optimx.run, optimx version 2014.5.4
  ## downloaded from r-forge

  # from top of optimx.run:
  # And make sure that controls NOT in method are deleted (nulled)
  mcontrol$follow.on <- NULL
  mcontrol$usenumDeriv <- NULL # JN130207
  mcontrol$save.failures <- NULL
  mcontrol$kkt <- NULL
  mcontrol$starttests <- NULL
  mcontrol$all.methods <- NULL
  mcontrol$dowarn <- NULL
  mcontrol$kkttol <- NULL
  mcontrol$kkt2tol <- NULL
  mcontrol$maximize <- NULL # Even if method DOES have it
  mcontrol$badval <- NULL
  mcontrol$scaletol <- NULL


  # from optimx.run nlminb section:
  # different name for iteration limit in this routine
  mcontrol$iter.max <- mcontrol$maxit
  mcontrol$maxit <- NULL
  mcontrol$abs.tol <- 0 # To fix issues when minimum is less than 0. 20100711
  if ( is.null(mcontrol$trace) || is.na(mcontrol$trace) || mcontrol$trace == 0){
    mcontrol$trace <- 0
  }else{
    # this is EVERY iteration. nlminb trace is freq of reporting.
    mcontrol$trace <- 1
  }

  # my addition here to get the parscaling
  if(!is.null(mcontrol$parscale)){
    scale <- 1/mcontrol$parscale
    mcontrol$parscale <- NULL
  }else{
    scale <- 1# NULL
  }
  time <- system.time(ans <- try(nlminb(start=par, objective=optcfg$ufn,
                                        lower=lower, upper=upper,
                                        control=mcontrol, scale=scale,
                                        ddfobj=ddfobj, ...),
                                 silent=TRUE))[1]
  if(class(ans)[1] != "try-error"){
    ans$convcode <- ans$convergence
    # Translate output to common format and names
    ans$value <- ans$objective
    ans$objective <- NULL
    ans$fevals <- ans$evaluations[1]
    ans$gevals <- ans$evaluations[2]
    ans$evaluations <- NULL # cleanup
    ans$nitns <- ans$iterations
    ans$iterations <- NULL
  }else{ # bad result -- What to do?
    ans <- list(fevals=NA) # ans not yet defined, so set as list
    ans$convcode <- 9999 # failed in run
    #if (ctrl$trace>0) cat("nlminb function evaluation failure\n")
    ans$value <- NA#ctrl$badval
    ans$par <- rep(NA,length(par))
    ans$nitns <- NA # not used
    ans$gevals <- NA ## ?? missing 130929
  }
  if(!exists("ans$message")) ans$message <- "none"
  ans$convergence <- NULL
  ans$xtimes <- time

  names(ans$par) <- names(par)

  cnames <- c(names(par), "value", "fevals", "gevals", "niter", "convcode",
              "kkt1", "kkt2", "xtimes")
  ans.ret <- matrix(NA, nrow=1, ncol=length(par)+8)
  colnames(ans.ret) <- cnames
  row.names(ans.ret) <- "nlminb"
  ans.ret["nlminb", ] <- c(ans$par, ans$value, ans$fevals, ans$gevals,
                           ans$nitns, ans$convcode, NA, NA, ans$xtimes)

  # from optimx.R again
  ans.details <- data.frame(method="nlminb", ngatend=NA, nhatend=NA,
                      hev=NA, message=ans$message)

  ansout <- data.frame(ans.ret)
  attr(ansout, "details") <- ans.details

  rownames(ans.details) <- "nlminb"
  # Fix kkt test output to logical
  ansout[ , "kkt1"] <- NA
  ansout[ , "kkt2"] <- NA

  # make an optimx object, so we can use optimx methods if we want
  structure(ansout, details = ans.details, maximize = FALSE,
            npar = length(par), follow.on=FALSE,
            class = c("optimx", "data.frame"))

}
