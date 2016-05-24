#' Plot a fitted \code{carx} object
#'
#' \code{plot.carx} plots a fitted \code{carx} object.
#' 
#' The y axis will be the values related to the response.
#' If the fitted object contains the data and censored  information in a \code{cenTS} object,
#'  the function will take advantage of the plot function for a \code{cenTS} object and superimpose the plot of the fitted values.
#'  Otherwise, the plot function will try to produce a plot similar to the previous case, while the x axis can be supplied by the user through \code{xAxisVar}, which must be ordinal and increasing.
#'
#' @param x a fitted \code{carx} object.
#' @param FUN an optional function to be applied to the transform the responses before plotting. 
#' This is useful for plotting the data on the original scale if the fitted \code{carx} object is based on
#' transformed responses. For instance, if the \code{carx} object was fitted with log transformed responses, 
#' setting \code{FUN} to exp renders the original response data to be plotted. Default = \code{NULL}.
#' @param xAxisVar an optional vector to be plotted as the x variable. Default =
#' \code{NULL} corresponds to doing a time series plot.
#' @param xlab the label of the x axis. Default = "Index".
#' @param ylab the label of the y axis. Default = "Response".
#' @param ... other parameters supplied to the generic function \code{plot}.
#' @return None. A plot will be displayed.
#'
#' @export
#' @examples
#'
#' #case 1: plot with cenTS object in the object, note that the x-axis is in date.
#' dat = carxSimCenTS(nObs=100,seed=0)
#' mdl <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)
#' #use default settings
#' plot(mdl)
#'
#' #case 2: plot without cenTS object in the object, note that the x-axis is a vector of numbers.
#' dat = carxSim(nObs=100,seed=0)
#' mdl <- carx(y=dat$y, x=dat[,c("X1","X2")], ci=dat$ci, lcl=dat$lcl, ucl=dat$ucl, p=2)
#' #or simply call
#' mdl <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)
#' plot(mdl)
plot.carx <- function(x,FUN=identity, xAxisVar=NULL, xlab="Index", ylab="Response",...)
{
	object <- x
  finiteRows <- object$finiteRows
  cts <- object$cenTS
  if(!is.null(cts) & is.null(xAxisVar))
  {
    cts$fitted = fitted(object)
    cts = FUN(cts)
    graphics::plot(cts,xlab=xlab,ylab=ylab,ylim=range(cts$fitted,na.rm=TRUE),...)
    graphics::lines(xts::.index(cts),cts$fitted,lty=2, col='blue')
    if(!is.null(object$outlier.indices))
      graphics::abline(v=xts::.index(cts)[object$outlier.indices],col="red",lty=2)
  }
  else
  {
    if(is.null(xAxisVar)) xAxisVar <- 1:length(object$y)
    ordinalX <- FALSE
    if(all(diff(xAxisVar)>0))
      ordinalX <- TRUE
    stopifnot(ordinalX)

    xrange <- range(xAxisVar)
    lgd <- c(ylab,'Fitted value')
    plty <-c(1,2)
    pcol <-c('black','blue')

    yh <- sapply(fitted(object), FUN)
    y <- sapply(object$y,  FUN)

    validLcl <- ifelse(any(is.finite(object$lcl)), TRUE, FALSE)
    validUcl <- ifelse(any(is.finite(object$ucl)), TRUE, FALSE)

    if(validLcl & ordinalX)
    {
      lcl <- sapply(object$lcl, FUN)
      lgd <- c(lgd,"Lower censoring limit")
      plty <- c(plty,3)
      pcol <- c(pcol,'red')
    }
    else
      lcl <- -Inf

    if(validUcl & ordinalX)
    {
      ucl <- sapply(object$ucl, FUN)
      lgd <- c(lgd,"Upper censoring limit")
      plty <- c(plty,5)
      pcol <- c(pcol,'red')
    }
    else
      ucl <- Inf

    ylim <- range(c(y,yh,lcl,ucl),na.rm=TRUE,finite=T)

    if( any(object$ci[finiteRows]>0) )
      y[finiteRows][object$ci[finiteRows]>0] <- ucl[finiteRows][object$ci[finiteRows]>0]
    if( any(object$ci[finiteRows]<0) )
      y[finiteRows][object$ci[finiteRows]<0] <- lcl[finiteRows][object$ci[finiteRows]<0]

    ylim[2] <- 1.3*ylim[2]

    graphics::plot(xAxisVar, y, type="l", lty=1, xlab=xlab, ylab=ylab, ylim=ylim, col='black',...)
    graphics::lines(xAxisVar, yh, lty=2, col='blue')

    if(validLcl & ordinalX)
    {
      graphics::lines(xAxisVar, lcl, lty=3, col="red")
      graphics::points(xAxisVar[object$ci<0],lcl[object$ci<0],pch=2)
    }
    if(validUcl & ordinalX)
    {
      graphics::lines(xAxisVar, ucl, lty=5, col="red")
      graphics::points(xAxisVar[object$ci>0], ucl[object$ci>0],pch=6)
    }
    graphics::legend("topright",legend=lgd,lty=plty,col=pcol)

    if(!is.null(object$outlier.indices))
      graphics::abline(v=xAxisVar[object$outlier.indices],col="red",lty=2)

  }
}

#' Show diagnostic plots for a \code{carx} object
#'
#' Four diagnostic plots will be shown, which are
#' \itemize{
#' \item{the time series plot of the residuals,}
#' \item{the residuals versus the fitted values,}
#' \item{the ACF plot of the residuals, and }
#' \item{the Ljung-Box test statistics versus the lags.}
#' }
#' @param object a \code{carx} object.
#' @param gof.lag  the maximum number of lags in ACF and Ljung-Box goodness-of-fit test.
#' @param col color of some warning lines in the figures, default=\code{red}.
#' @param omit.initial whether initial residuals should be omitted, default = \code{TRUE}.
#' @param mfrow \code{par} parameter indicating how the plots are to be arranged, default = \code{c(4,1)}.
#' @param main The main title of the plot, default = "Diagnostic Plots".
#' @param ... Other arguments sent to \code{plot}.
#' @return none.
#' @export
#' @examples
#' dat = carxSim(nObs=100,seed=0)
#' mdl <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)
#' tsdiag(mdl)


tsdiag.carx <- function(object
                        ,gof.lag
                        ,col="red"
                        ,omit.initial=TRUE
                        ,mfrow=c(4,1)
                        ,main="Diagnostic Plots"
                        #,plotRes=TRUE
                        #,plotResVSFit=TRUE
                        #,plotACF=TRUE
                        #,plotLB=TRUE
                        ,...)
{
  #nPlot <- sum(plotRes + plotResVSFit + plotACF + plotLB)

  opar = graphics::par(mfrow = mfrow, mar = c(4, 4, 4, 3) + 0.1, oma = c(1, 0, 2, 0))
  n = object$nObs
  cts = object$cenTS
  if (missing(gof.lag))
      lag.max = 10 * log10(n)
  else
    lag.max = gof.lag

  residuals = residuals(object,type="raw")
  f = fitted(object)
  if(omit.initial)
  {
    residuals = window(residuals, start = stats::time(residuals)[object$p + 1])
    f = window(f, start = stats::time(f)[object$p + 1])
    if(!is.null(cts))
      cts = window(cts, start = stats::time(cts)[object$p + 1])
  }

  std.res = residuals/object$sigma
  n = length(std.res)
  h1 = stats::qnorm(0.025/n)
  if(is.null(cts))
    graphics::plot(std.res, ylab = "Standardized Residuals", type = "p",main=main, ...)
  else
    graphics::plot(xts::xts(std.res,zoo::index(cts)),ylab = "Standardized Residuals", type = "p", main=main,...)
  graphics::abline(h = h1, lty = 2, col = col)
  graphics::abline(h = -h1, lty = 2, col = col)
  graphics::abline(h = 0)


  graphics::plot(f,std.res,xlab="Fitted",ylab="Standardized Residuals")

  stats::acf(as.numeric(std.res), lag.max = lag.max, ylab = "ACF of Residuals",
      ci.col = col, main = "", na.action=stats::na.pass,...)
  lbv = rep(NA, lag.max)
  k = object$p+1
  for (i in k:lag.max) {
      lbv[i] = stats::Box.test(std.res, lag = i, type="Ljung-Box",fitdf = object$p)$p.value
  }
  graphics::plot(y = lbv, x = 1:lag.max, ylim = c(0, 1), pch = 21, ylab = "P-values",
      xlab = "Number of lags", ...)
  graphics::abline(h = 0.05, lty = 2, col = col)
  graphics::par(opar)
  invisible()
}
