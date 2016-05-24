#' Randomised quantile residuals check plots for GAMs/DSMs
#'
#' Function operates as \code{\link{gam.check}} but using randomised quantile
#' residuals, a la Dunn and Smyth (1996). Checks of \code{k} are not computed,
#' these need to be done using \code{\link{gam.check}}.
#'
#' @param gam.obj a \code{gam}, \code{glm} or \code{dsm} object.
#' @param ... arguments passed on to all plotting functions
#' @return just plots!
#'
#' @author Based on code provided by Natalie Kelly, bugs added by Dave Miller
#' @export
#' @importFrom statmod qres.nbinom qresid qres.tweedie
#' @importFrom graphics par hist
#' @importFrom stats napredict fitted qqnorm
#'
#' @examples
#' \donttest{
#' library(Distance)
#' library(dsm)
#' library(tweedie)
#'
#' # load the Gulf of Mexico dolphin data (see ?mexdolphins)
#' data(mexdolphins)
#'
#' # fit a detection function and look at the summary
#' hr.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance),
#'                key = "hr", adjustment = NULL)
#'
#' # fit a simple smooth of x and y with a Tweedie response with estimated
#' #  p parameter
#' mod1<-dsm(N~s(x,y), hr.model, mexdolphins$segdata, mexdolphins$obsdata, family=tw())
#' rqgam.check(mod1)
#' }
rqgam.check<-function(gam.obj,...){

  # layout stuff
  opar <- par(mfrow=c(2,2))

  # grab the randomised quantile residuals
  # requires statmod package

  # need to do the right thing for mgcv's Tweedie
  if(grepl("^Tweedie",gam.obj$family$family)){
    if(is.null(environment(gam.obj$family$variance)$p)){
      p.val <- gam.obj$family$getTheta(TRUE)
      environment(gam.obj$family$variance)$p <- p.val
    }
    qres <- qres.tweedie(gam.obj)
  # and for negbin
  }else if(grepl("^Negative Binomial",gam.obj$family$family)){
    # need to set $theta
    if("extended.family" %in% class(gam.obj$family)){
      # for SNW's extended family, need to set TRUE in getTheta as theta
      # is on the wrong scale
      gam.obj$theta <- gam.obj$family$getTheta(TRUE)
    }else{
      gam.obj$theta <- gam.obj$family$getTheta()
    }
    qres <- qres.nbinom(gam.obj)
  }else{
    # for everything else
    qres <- qresid(gam.obj)
  }

  # values of the linear predictor
  linpred <- napredict(gam.obj$na.action, gam.obj$linear.predictors)

  ## normal Q-Q plot
  qqnorm(qres,ylab="Randomised quantile residuals",...)

  ## resids vs. linear pred
  plot(linpred, qres,main="Resids vs. linear pred.",
         xlab="linear predictor",ylab="Randomised quantile residuals",...)

  ## histogram
  hist(qres, main="Histogram of residuals",
       xlab="Randomised quantile residuals",...)


  ## Response vs. Fitted Values
  plot(fitted(gam.obj), gam.obj$y,
       main="Response vs. Fitted Values",
       xlab="Fitted Values", ylab="Response",...)
  #lines(lowess(gam.obj$fitted.values, gam.obj$model[,1]), col = 2)


  #residuals versus leverage
  #plot(gam.obj$hat, qresid(gam.obj), main="QRES versus leverage", xlab="Leverage", ylab="QRES")
  #lines(lowess(gam.obj$hat, qresid(gam.obj)), col = 2)
  #plot(gam.obj$fitted.values, sqrt(abs(qres)),
  #     main="Scale-Location (QRES)",
  #     xlab="Fitted", ylab="Sqrt(Abs(QRES))" )
  ##lines(lowess(gam.obj$fitted.values, sqrt(abs(qresid(gam.obj)))), col = 2)


  par(opar)
}
