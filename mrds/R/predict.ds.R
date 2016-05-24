#' Predictions from \code{mrds} models
#'
#' Predict detection probabilities (or effective strip widths/effective areas of detection) from a fitted distance sampling model using either the original data (i.e. "fitted" values) or using new data.
#'
#' The first 4 arguments are the same in each predict function.  The latter 2 are specific to certain functions. For line transects, the effective strip half-width (\code{esw=TRUE}) is the integral of the fitted detection function over either 0 to W or the specified \code{int.range}.  The predicted detection probability is the average probability which is simply the integral divided by the distance range. For point transect models, \code{esw=TRUE} calculates the effective area of detection (commonly referred to as "nu", this is the integral of \code{2/width^2 * rg(r)}.
#'
#' Fitted detection probabilities are stored in the \code{model} object and these are returned unless \code{compute=TRUE} or \code{newdata} is specified. \code{compute=TRUE} is used to estimate numerical derivatives for use in delta method approximations to the variance.
#'
#' For \code{method="io.fi"} or \code{method="trial.fi"} if
#' \code{integrate=FALSE}, \code{predict} returns the value of the conditional
#'  detection probability and if \code{integrate=TRUE}, it returns the average
#' conditional detection probability by integrating over x (distance) with
#' respect to a uniform distribution.
#'
#' @aliases predict predict.ds predict.ddf predict.io predict.io.fi predict.trial
#'  predict.trial.fi predict.rem predict.rem.fi
#' @param object \code{ddf} model object.
#' @param newdata new \code{data.frame} for prediction.
#' @param compute if \code{TRUE} compute values and don't use the fitted
#'  values stored in the model object.
#' @param int.range integration range for variable range analysis; either
#'  vector or matrix.
#' @param esw if \code{TRUE}, returns effective strip half-width (or effective area of detection for point transect models) integral from 0 to the truncation distance (\code{width}) of p(y)dy; otherwise it returns the integral from 0 to truncation width of p(y)*pi(y) where pi(y)=1/W for lines and pi(y)=2r/W^2 for points.
#' @param integrate for \code{*.fi} methods, see Details below.
#' @param \dots for S3 consistency
#' @usage \method{predict}{ds}(object,newdata,compute=FALSE,int.range=NULL,esw=FALSE,...)
#'        \method{predict}{io.fi}(object,newdata,compute=FALSE, int.range=NULL,integrate=FALSE,...)
#'        \method{predict}{io}(object,newdata,compute=FALSE,int.range=NULL,...)
#'        \method{predict}{trial}(object,newdata,compute=FALSE,int.range=NULL,...)
#'        \method{predict}{trial.fi}(object,newdata,compute=FALSE, int.range=NULL,integrate=FALSE,...)
#'        \method{predict}{rem}(object,newdata,compute=FALSE,int.range=NULL,...)
#'        \method{predict}{rem.fi}(object,newdata,compute=FALSE, int.range=NULL,integrate=FALSE,...)
#' @return For all but the exceptions below, the value is a list with a single
#'   element: \code{fitted}, a vector of average detection probabilities or esw values for each observation in the original data or\code{newdata}
#'
#' For \code{predict.io.fi},\code{predict.trial.fi},\code{predict.rem.fi} with \code{integrate=TRUE}, the value is a list with one element: \code{fitted}, which is a vector of integrated (average) detection probabilities for each observation in the original data or \code{newdata}.
#'
#' For \code{predict.io.fi}, \code{predict.trial.fi}, or \code{predict.rem.fi}
#'   with \code{integrate=FALSE}, the value is a list with the following
#'   elements:
#'  \describe{
#'    \item{\code{fitted}}{\eqn{p(y)} values}
#'    \item{\code{p1}}{\eqn{p_{1|2}(y)}, conditional detection probability for observer 1}
#'    \item{\code{p2}}{\eqn{p_{2|1}(y)}, conditional detection probability for observer 2}
#'    \item{\code{fitted}}{\eqn{p_.(y)=p_{1|2}(y)+p_{2|1}(y)-p_{1|2}(y)*p_{2|1}(y)}, conditional detection probability of being seen by either observer}}
#'
#' @note Each function is called by the generic function \code{predict} for the
#'   appropriate \code{ddf} model object.  They can be called directly by the
#'   user, but it is typically safest to use \code{predict} which calls the
#'   appropriate function based on the type of model.
#' @author Jeff Laake
#' @seealso \code{\link{ddf}}, \code{\link{summary.ds}},
#'   \code{\link{plot.ds}}
#' @keywords utility
#' @export
#' @importFrom stats predict
# Uses: integratedetfct, integratedetfct.logistic
predict.ds <- function(object,newdata=NULL,compute=FALSE,int.range=NULL,
                       esw=FALSE,...){
  model <- object
  ltmodel <- model$ds
  x <- ltmodel$aux$ddfobj$xmat
  point <- ltmodel$aux$point
  width <- ltmodel$aux$width

  # If there are no fitted values present or compute is set TRUE or
  # a newdata frame has been used, then the predicted values must be computed.
  if(is.null(model$fitted) | compute | !is.null(newdata)){

    # Get model and extract the parameters. Note that in computing derivatives
    # for variances, the model parameters are perturbed and may not be at
    # mle values after fitting.
    fpar <- model$par
    ddfobj <- ltmodel$aux$ddfobj
    ddfobj <- assign.par(ddfobj,fpar)

    # Get integration ranges either from specified argument or from
    # values stored in the model.
    if(is.null(int.range)){
      if(is.null(newdata)){
        nr <- nrow(ddfobj$xmat)
      }else{
        nr <- nrow(newdata)
      }

      if(is.null(ltmodel$aux$int.range)){
        int.range <- cbind(rep(0,nr),rep(width,nr))
      }else{
        int.range <- ltmodel$aux$int.range
        if(is.vector(int.range)){
          int.range <- cbind(rep(int.range[1],nr),
                             rep(int.range[2],nr))
        }
      }
    }

    # Extract other values from model object
    if(!is.null(newdata)){
      if(!is.null(ddfobj$scale)){
        zdim <- ncol(ddfobj$scale$dm)
        znames <- colnames(ddfobj$scale$dm)
        ddfobj$scale$dm <- setcov(newdata, as.formula(ddfobj$scale$formula))
        if(zdim != ncol(ddfobj$scale$dm) |
           !all(znames==colnames(ddfobj$scale$dm)) ){
          stop("fields or factor levels in newdata do not match data used in estimation model for scale model\n")
        }
      }

      if(!is.null(ddfobj$shape)){
        zdim <- ncol(ddfobj$shape$dm)
        znames <- colnames(ddfobj$shape$dm)
        ddfobj$shape$dm <- setcov(newdata, as.formula(ddfobj$shape$formula))
        if(zdim != ncol(ddfobj$shape$dm) |
           !all(znames==colnames(ddfobj$shape$dm))){
          stop("fields or factor levels in newdata do not match data used in estimation model for shape model\n")
        }
      }
      # update xmat too
      datalist <- process.data(newdata,object$meta.data,check=FALSE)
      ddfobj$xmat <- datalist$xmat
    }

    # Compute integral of fitted detection function using either logistic or
    # non-logistic detection function.  Note that "logistic" is not currently
    # allowed as it has not been fully tested.

    # if(ftype=="logistic")
    #   int1=integratedetfct.logistic(x,ltmodel$model$scalemodel,width,
    #                         int.range,theta1,ltmodel$aux$integral.numeric,z)
    # else
    int1 <- integratepdf(ddfobj,select=rep(TRUE,nrow(ddfobj$xmat)),width=width,
                         int.range=int.range,
                         standardize=TRUE,point=point)
    # int1=integratedetfct(ddfobj,select=rep(TRUE,nrow(ddfobj$xmat)),
    #           width=width,int.range=int.range,point=point)

    # If the predicted values don't need to be computed, then use the values
    # in the model object (model$fitted) and change to integral (esw) values.
    # Note this needs to be checked to see if it works with variable ranges.

  }else{
    int1 <- model$fitted
  }

  # Compute either esw (int1) or p and store in fitted.
  if(esw){
    if(!point){
      fitted <- int1*width
    }else{
      fitted <- int1*pi*width^2
    }
  }else{
     fitted <- int1
  }

  # If there are no covariates and there is only one prediction, expand to
  # a vector based on length of data object.
  if(length(fitted)==1){
    if(is.null(newdata)){
      fitted <- rep(fitted,length(x$object))
    }else{
      fitted <- rep(fitted,nrow(newdata))
    }
  }

  # If not a new dataframe, then assign names from data stored in the model
  # object otherwise, use those from newdata.  Then return vector of values
  # to calling frame.
  if(is.null(newdata)){
    names(fitted) <- x$object
  }else{
    names(fitted) <- newdata$object
  }

  return(list(fitted=fitted))
}
