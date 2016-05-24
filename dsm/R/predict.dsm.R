#' Predict from a fitted density surface model
#'
#' Make predictions outside (or inside) the covered area.
#'
#' @param object a fitted \code{\link{dsm}} object as produced by \code{dsm()}.
#' @param newdata spatially referenced covariates e.g. altitude, depth, distance to shore, etc. Covariates in the \code{data.frame} must have names *identical* to variable names used in fitting the DSM.
#' @param off.set area of each of the cells in the prediction grid. Should be in the same units as the segments/distances given to \code{dsm}. Ignored if there is already a column in \code{newdata} called \code{off.set}.
#' @param type what scale should the results be on. The default is
#'  \code{"response"}, see \code{\link{predict.gam}} for an explanation of other options (usually not necessary).
#' @param \dots any other arguments passed to \code{\link{predict.gam}}.
#' @return predicted values on the response scale (density/abundance).
#' @export
#'
#' @seealso predict.gam dsm.var.gam dsm.var.prop dsm.var.movblk
#' @author David L. Miller
#' @importFrom stats predict
predict.dsm <- function(object, newdata=NULL, off.set=NULL,
                        type="response",...){

  if("gamm" %in% class(object)){
    object <- object$gam
  }

  if(is.null(newdata)){
    newdata <- object$data
  }

  # if we don't have a density model, then set the offset
  if(!(c(object$formula[[2]]) %in% c("D","presence","density"))){
    if(is.null(newdata$off.set)){
      if(is.null(off.set)){
        stop("You must supply off.set in data or as an argument.")
      }
      newdata$off.set <- off.set
    }

    # apply the link function
    linkfn <- object$family$linkfun
    newdata$off.set <- linkfn(newdata$off.set)
  }

  # remove the dsm class
  class(object) <- class(object)[class(object)!="dsm"]

  # actually do the predict call
  result<-predict(object, newdata, type=type,...)

  return(result)
}
