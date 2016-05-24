#' Variance estimation via Bayesian results
#'
#' Use results from the Bayesian interpretation of the GAM to obtain uncertainty estimates. See Wood (2006).
#'
#' This is based on \code{\link{dsm.var.prop}} taken from code by Mark Bravington and Sharon Hedley.
#'
#' @param dsm.obj an object returned from running \code{\link{dsm}}.
#' @param pred.data either: a single prediction grid or list of prediction grids. Each grid should be a \code{data.frame} with the same columns as the original data.
#' @param off.set a a vector or list of vectors with as many elements as there are in \code{pred.data}. Each vector is as long as the number of rows in the corresponding element of \code{pred.data}. These give the area associated with each prediction cell. If a single number is supplied it will be replicated for the length of \code{pred.data}.
#' @param seglen.varname name for the column which holds the segment length (default value \code{"Effort"}).
#' @param type.pred should the predictions be on the "response" or "link" scale? (default \code{"response"}).
#' @return a list with elements
#'         \tabular{ll}{\code{model} \tab the fitted model object\cr
#'                      \code{pred.var} \tab variance of the regions given
#'                      in \code{pred.data}.\cr
#'                      \code{bootstrap} \tab logical, always \code{FALSE}\cr
#'                      \code{model}\tab the fitted model with the extra term\cr
#'                      \code{dsm.object} \tab the original model, as above
#'                      }
#' @author David L. Miller
#' @importFrom stats coef vcov
#' @export
dsm.var.gam<-function(dsm.obj, pred.data,off.set=NULL,
    seglen.varname='Effort', type.pred="response") {

  # strip dsm class so we can use gam methods
  class(dsm.obj) <- class(dsm.obj)[class(dsm.obj)!="dsm"]

  # if we have a gamm, then just pull out the gam object
  if(any(class(dsm.obj)=="gamm")){
    dsm.obj <- dsm.obj$gam
    is.gamm <- TRUE
  }

  # if all the offsets are the same then we can just supply 1 and rep it
  if(length(off.set)==1){
    if(is.null(nrow(pred.data))){
      off.set <- rep(list(off.set),length(pred.data))
    }else{
      off.set <- rep(off.set,nrow(pred.data))
    }
  }

  # make sure if one of pred.data and off.set is not a list we break
  # if we didn't have a list, then put them in a list so everything works
  if(is.data.frame(pred.data) & is.vector(off.set)){
    pred.data <- list(pred.data)
    off.set <- list(off.set)
#    pred.data[[1]] <- pred.data
#    off.set[[1]] <- off.set
  }else if(is.list(off.set)){
    if(length(pred.data)!=length(off.set)){
      stop("pred.data and off.set don't have the same number of elements")
    }
  }

  # depending on whether we have response or link scale predictions...
  if(type.pred=="response"){
      tmfn <- dsm.obj$family$linkinv
      dtmfn <- function(eta){vapply(eta, numderiv, numeric(1), f=tmfn)}
  }else if(type.pred=="link"){
      tmfn <- identity
      dtmfn <- function(eta){1}
  }

  # grab the coefficients
  cft <- coef(dsm.obj)
  preddo <- list(length(pred.data))
  dpred.db <- matrix(0, length(pred.data), length(cft))

  # loop over the prediction grids
  for(ipg in seq_along(pred.data)){
    ### fancy lp matrix stuff
    # set the offset to be zero here so we can use lp
    pred.data[[ipg]]$off.set<-rep(0,nrow(pred.data[[ipg]]))

    lpmat <- predict(dsm.obj, newdata=pred.data[[ ipg]], type='lpmatrix')
    lppred <- lpmat %**% cft

    # if the offset is just one number then repeat it enough times
    if(length(off.set[[ipg]])==1){
      this.off.set <- rep(off.set[[ipg]],nrow(pred.data[[ipg]]))
    }else{
      this.off.set <- off.set[[ipg]]
    }

    preddo[[ipg]] <-  this.off.set %**% tmfn(lppred)
    dpred.db[ipg,] <- this.off.set %**% (dtmfn(lppred)*lpmat)
  }

  # "'vpred' is the covariance of all the summary-things." - MVB
  # so we want the diagonals if length(pred.data)>1
  # A B A^tr
  vpred <- dpred.db %**% tcrossprod(vcov(dsm.obj), dpred.db)

  if(is.matrix(vpred)){
    vpred <- diag(vpred)
  }

  result <- list(pred.var       = vpred,
                 bootstrap      = FALSE,
                 pred           = preddo,
                 var.prop       = FALSE,
                 pred.data      = pred.data,
                 off.set        = off.set,
                 dsm.object     = dsm.obj,
                 seglen.varname = seglen.varname,
                 type.pred      = type.pred
                )

  class(result) <- "dsm.var"

  return(result)
}
