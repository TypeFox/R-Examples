#' Extract coefficients
#'
#' Extract coefficients and provide a summary of parameters and estimates
#' from the output of \code{\link{ddf}} model objects.
#'
#' @aliases coefficients coef.ds coef.io coef.io.fi coef.trial coef.trial.fi
#'   coef.rem coef.rem.fi
#' @param object ddf model object of class \code{ds}, \code{io}, \code{io.fi},
#'    \code{trial}, \code{trial.fi}, \code{rem}, or \code{rem.fi}.
#' @param \dots unspecified arguments that are unused at present
#' @usage \method{coef}{ds}(object,...)
#'        \method{coef}{io}(object,...)
#'        \method{coef}{io.fi}(object,...)
#'        \method{coef}{trial}(object,...)
#'        \method{coef}{trial.fi}(object,...)
#'        \method{coef}{rem}(object,...)
#'        \method{coef}{rem.fi}(object,...)
#' @export
#' @return
#'
#' For \code{coef.ds} List of data frames for coefficients (scale and exponent
#'   (if hazard)) \item{scale}{dataframe of scale coefficent estimates and
#'   standard errors} \item{exponent}{dataframe with exponent estimate and
#'   standard error if hazard detection function}
#'
#' For all others Data frame containing each coefficient and standard error
#' @note These functions are called by the generic function \code{coef} for any
#'   \code{ddf} model object.  It can be called directly by the user, but it is
#'   typically safest to use \code{coef} which calls the appropriate function
#'   based on the type of model.
#' @author Jeff Laake
#' @importFrom stats coef
coef.ds <- function(object,...){
  ltmodel <- object$ds
  coeff <- NULL
  if(is.null(object$par)){
    vcov <- NULL
  }else{
    vcov <- solvecov(object$hessian)$inv
  }

  ddfobj <- ltmodel$aux$ddfobj
  indices <- getpar(ddfobj,index=TRUE)
  se <- sqrt(diag(vcov))

  if(indices[1]!=0){
    key.shape.se <- se[1:indices[1]]
  }

  if(indices[2]!=0){
    key.scale.se <- se[(indices[1]+1):indices[2]]
  }

  if(indices[3]!=0){
    if(indices[2]!=0){
      adj.parm.se <- se[(indices[2]+1):indices[3]]
    }else{
      adj.parm.se <- se[(indices[1]+1):indices[3]]
    }
  }

  # Always get the scale parameter
  if(!is.null(ddfobj$scale)){
    coeff <- as.data.frame(cbind(estimate=ddfobj$scale$parameters,
                                 se=key.scale.se))
    row.names(coeff) <- colnames(ddfobj$scale$dm)
  }

  # If we have a non-null shape parameter, get that too
  if(!is.null(ddfobj$shape)){
    exp.coeff <- as.data.frame(cbind(estimate=ddfobj$shape$parameters,
                                     se=key.shape.se))
    row.names(exp.coeff) <- colnames(ddfobj$shape$dm)
  }

  # Get the adjustment parameters if necessary
  if(!is.null(ddfobj$adjustment)){
    adj.coeff <- as.data.frame(cbind(estimate=ddfobj$adjustment$parameters,
                                     se=adj.parm.se))
    adj.names <- NULL
    for(i in 1:nrow(adj.coeff)){
      adj.names[i] <- paste(ddfobj$adjustment$series,
                            ", order ",ddfobj$adjustment$order[i],sep="")
    }
    row.names(adj.coeff) <- adj.names

    # Return values if we have adjustment terms
    if(!is.null(ddfobj$shape)){
      return(list(scale=coeff,exponent=exp.coeff,adjustment=adj.coeff))
    }else if(is.null(ddfobj$scale)){
      return(list(adjustment=adj.coeff))
    }else{
      return(list(scale=coeff,adjustment=adj.coeff))
    }

  }

  # Return values if no adjustment terms
  if(!is.null(ddfobj$shape)){
    return(list(scale=coeff,exponent=exp.coeff))
  }else{
    return(list(scale=coeff))
  }
}
