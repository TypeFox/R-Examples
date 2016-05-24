#' Estimate parameters from a specified model using estimated maximum likelihood
#'
#' @param start Vector of starting values, if NULL, will come up with starting
#'   values
#' @param method Method to use for optimization, can be "pseudo-score" for
#'   categorical BIP, or any of the methods
#'   available in \link{optim}. Defaults to "BFGS"
#' @param control List of control parameters for passed to \link{optim}
#' @param ... Arguments passed to \link{optim} or \link{pseudo_score}.
#'
#' @export
#'
ps_estimate <- function(start = NULL, method = "BFGS", control = list(), ...){

  est.call <- match.call()
  rval <- function(psdesign){

    if(!"risk.model" %in% names(psdesign)) stop("No risk model specified")
    if(!"integration.models" %in% names(psdesign)) stop("No integration models specified")


    if(is.null(start)){

      start <- rep(0, psdesign$nparam)

    }

    if(method == "pseudo-score"){

      est1 <- pseudo_score(psdesign, start = start, ...)
      names(est1$par) <- psdesign$param.names


    } else {

      est1 <- optim(start, fn = psdesign$likelihood, method = method, control = control, ...)
      names(est1$par) <- psdesign$param.names

    }
    psdesign$estimate.call <- est.call
    psdesign$estimates <- est1
    psdesign

  }

  class(rval) <- c("ps", "estimate")
  rval

}


#' Estimate parameters from a specified model using bootstrap resampling and estimated maximum likelihood
#'
#' @param n.boots Number of bootstrap replicates
#' @param progress.bar Logical, if true will display a progress bar in the console
#' @param start Vector of starting values, if NULL, will come up with starting values
#' @param method Method to use for optimization, can be "pseudo-score" for
#'   categorical S with nonparametric integration, or any of the methods
#'   available in \link{optim}. Defaults to "BFGS"
#' @param control List of control parameters for passed to \link{optim}
#' @param ... Arguments passed to \link{optim}
#'
#' @export
#


ps_bootstrap <- function(n.boots = 200, progress.bar = TRUE, start = NULL, method = "BFGS", control = list(), ...){

  bs.call <- match.call()

  rval <- function(psdesign){
    if(!"risk.model" %in% names(psdesign)) stop("No risk model specified")
    if(!"integration.models" %in% names(psdesign)) stop("No integration models specified")

    bootpar <- vector(mode = "list", length = n.boots)

    if(progress.bar){
      cat(paste("Bootstrapping", n.boots, "replicates:\n"))
      pb <- txtProgressBar(min = 1, max = n.boots)
    }

    psdesign.0 <- psdesign
    for(i in 1:n.boots){
      # resample augdata stratified by Y


        Z0 <- (1:nrow(psdesign$augdata))[psdesign$augdata$Z == 0]
        Z1 <- (1:nrow(psdesign$augdata))[psdesign$augdata$Z == 1]


      sampdex <- c(sample(Z0, length(Z0), replace = TRUE), sample(Z1, length(Z1), replace = TRUE))
      psdesign.0$augdata <- psdesign$augdata[sampdex, ]

      if(is.factor(psdesign.0$augdata$S.1)){

        if(length(unique(as.numeric(psdesign.0$augdata$S.1))) != length(unique(as.numeric(psdesign$augdata$S.1)))){
          bootpar[[i]] <- c(rep(NA, psdesign$nparam), convergence = 11)
          next
        }

      }
      ## re-call integration models

      psdesign2 <- psdesign.0
      for(intj in psdesign$integration.models){
        psdesign2 <- psdesign2 + do.call(as.character(intj$model$args[[1]]), intj$model$args[-1])
      }

      ## re-call risk model

      psdesign3 <- psdesign2 + do.call(as.character(psdesign$risk.model$args[[1]]), psdesign$risk.model$args[-1])

      # estimate

      bpar <- psdesign3 + ps_estimate(start = start, control = control, ...)

      bootpar[[i]] <- c(bpar$estimates$par, convergence = bpar$estimates$convergence)

      if(progress.bar){
        setTxtProgressBar(pb, value = i)
        flush.console()
      }

    }

    if(progress.bar) close(pb)

    bootpar <- as.data.frame(do.call(rbind, bootpar))
    psdesign$bootstraps <- bootpar
    psdesign$bs.call <- bs.call
    psdesign

  }

  class(rval) <- c("ps", "bootstrap")
  rval

}

