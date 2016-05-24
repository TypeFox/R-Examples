#' Plotting predicted survival curves.
#' 
#' Ploting prediction survival curves for one prediction model using
#' \code{predictSurvProb} .
#' 
#' Arguments for the invoked functions \code{legend} and \code{axis} are simply
#' specified as \code{legend.lty=2}. The specification is not case sensitive,
#' thus \code{Legend.lty=2} or \code{LEGEND.lty=2} will have the same effect.
#' The function \code{axis} is called twice, and arguments of the form
#' \code{axis1.labels}, \code{axis1.at} are used for the time axis whereas
#' \code{axis2.pos}, \code{axis1.labels}, etc.  are used for the y-axis.
#' 
#' These arguments are processed via \code{\dots{}} of
#' \code{plotPredictSurvProb} and inside by using the function
#' \code{SmartControl}.
#' 
#' @param x A survival prediction model including \code{call} and
#' \code{formula} object.
#' @param newdata A data frame with the same variable names as those that were
#' used to fit the model \code{x}.
#' @param times Vector of times at which to return the estimated probabilities.
#' @param xlim Plotting range on the x-axis.
#' @param ylim Plotting range on the y-axis.
#' @param xlab Label given to the x-axis.
#' @param ylab Label given to the y-axis.
#' @param axes Logical. If \code{FALSE} no axes are drawn.
#' @param col Vector of colors given to the survival curve.
#' @param density Densitiy of the color -- useful for showing many
#' (overlapping) curves.
#' @param lty Vector of lty's given to the survival curve.
#' @param lwd Vector of lwd's given to the survival curve.
#' @param add Logical. If \code{TRUE} only lines are added to an existing
#' device
#' @param legend Logical. If TRUE a legend is plotted by calling the function
#' legend.  Optional arguments of the function \code{legend} can be given in
#' the form \code{legend.x=val} where x is the name of the argument and val the
#' desired value. See also Details.
#' @param percent Logical. If \code{TRUE} the y-axis is labeled in percent.
#' @param \dots Parameters that are filtered by \code{\link{SmartControl}} and
#' then passed to the functions: \code{\link{plot}}, \code{\link{axis}},
#' \code{\link{legend}}.
#' @return The (invisible) object.
#' @author Ulla B. Mogensen \email{ulmo@@biostat.ku.dk}, Thomas A. Gerds
#' \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{predictSurvProb}}\code{\link{prodlim}}
#' @references Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
#' Evaluating Random Forests for Survival Analysis Using Prediction Error
#' Curves. Journal of Statistical Software, 50(11), 1-23. URL
#' http://www.jstatsoft.org/v50/i11/.
#' @keywords survival
#' @examples
#' 
#' # generate some survival data
#' library(prodlim)
#' d <- SimSurv(100)
#' # then fit a Cox model 
#' library(rms)
#' coxmodel <- cph(Surv(time,status)~X1+X2,data=d,surv=TRUE)
#' # plot predicted survival probabilities for all time points
#' ttt <- sort(unique(d$time))
#' # and for selected predictor values:
#'  ndat <- data.frame(X1=c(0.25,0.25,-0.05,0.05),X2=c(0,1,0,1))
#' plotPredictSurvProb(coxmodel,newdata=ndat,times=ttt)
#' 
#' # the same can be done e.g. for a randomSurvivalForest model
#' library(randomForestSRC)
#' rsfmodel <- rfsrc(Surv(time,status)~X1+X2,data=d)
#' plotPredictSurvProb(rsfmodel,newdata=ndat,times=ttt)
#' 
#' @export
plotPredictSurvProb <- function(x,
                                newdata,
                                times,
                                xlim,
                                ylim,
                                xlab,
                                ylab,
                                axes=TRUE,
                                col,
                                density,
                                lty,
                                lwd,
                                add=FALSE,
                                legend=TRUE,
                                percent=FALSE,
                                ...){
  # {{{ call argument

  allArgs <- match.call()
  
  # }}}
  # {{{ find times

  if(missing(times)){
    # formula
    formula <- eval(x$call$formula)
    if (match("formula",class(formula),nomatch=0)==0)
      stop("Argument formula is missing.")
    # find data
    data <- eval(x$call$data)
    # extract response
    m <- model.frame(formula,data,na.action=na.fail)
    response <- model.response(m)
    # ordering time 
    neworder <- order(response[,"time"],-response[,"status"])
    response <- response[neworder,,drop=FALSE]
    times <- response[,"time"]
    # unique event times
    times <- unique(times)
  }

  # }}}
  # {{{ newdata
  if(missing(newdata)){
    newdata <- eval(x$call$data)
  }
  ## stop("newdata argument is missing")

  # }}}
  # {{{ xlim, ylim

  if (missing(xlim)) xlim <- c(0, max(times))
  at <- times <= xlim[2]
  orig.X <- times[at]
  X <- times[at]
  
  # }}}

  # {{{ predict newdata at times
  y <- predictSurvProb(x, newdata=newdata, times=orig.X)
   
  # }}}
  # {{{ plot arguments

  nlines <- NROW(y)
  
  if (missing(ylab)) ylab <- "Survival probability"
  if (missing(xlab)) xlab <- "Time"
  if (missing(ylim)) ylim <- c(0, 1)
  if (missing(lwd)) lwd <- rep(3,nlines)
  if (missing(col)) col <- rep(1,nlines)
  if (missing(density)){
    if (nlines>5){
      density <- pmax(33,100-nlines)
    }
    else density <- 100
  }
  if (density<100){
    col <- sapply(col,function(coli){
      ccrgb=as.list(col2rgb(coli,alpha=TRUE))
      names(ccrgb) <- c("red","green","blue","alpha")
      ccrgb$alpha=density
      cc=do.call("rgb",c(ccrgb,list(max=255)))
    })
  }
  if (missing(lty)) lty <- rep(1, nlines)
  if (length(lwd) < nlines) lwd <- rep(lwd, nlines)
  if (length(lty) < nlines) lty <- rep(lty, nlines)
  if (length(col) < nlines) col <- rep(col, nlines)
  
  axis1.DefaultArgs <- list()
  
  axis2.DefaultArgs <- list(at=seq(0,1,.25))
         
  plot.DefaultArgs <- list(x=0,
                           y=0,
                           type = "n",
                           ylim = ylim,
                           xlim = xlim,
                           xlab = xlab,
                           ylab = ylab)
  
  
  
  legend.DefaultArgs <- list(legend=rownames(y),
                             lwd=lwd,
                             col=col,
                             lty=lty,
                             cex=1.5,
                             bty="n",
                             y.intersp=1.3,
                             x="topright")
  
  # }}}
  # {{{ smart controls
 
  if (match("legend.args",names(args),nomatch=FALSE)){
    legend.DefaultArgs <- c(args[[match("legend.args",names(args),nomatch=FALSE)]],legend.DefaultArgs)
    legend.DefaultArgs <- legend.DefaultArgs[!duplicated(names(legend.DefaultArgs))]
  }
  smartA <- prodlim::SmartControl(call=list(...),
                                   keys=c("plot","legend","axis1","axis2"),
                                   ignore=c("x", "newdata", "times", "xlim","ylim","xlab","ylab","col","lty","lwd","add","legend","percent","axes","legend.args"),
                                   defaults=list("plot"=plot.DefaultArgs,
                                     "legend"= legend.DefaultArgs,
                                     "axis1"=axis1.DefaultArgs,
                                     "axis2"=axis2.DefaultArgs),
                                   forced=list("plot"=list(axes=FALSE),
                                     "axis1"=list(side=1),
                                     "axis2"=list(side=2)),
                                   verbose=TRUE)

  # }}} 
  # {{{ empty plot
  
  if (!add) {
    do.call("plot",smartA$plot)
    if (axes){
      do.call("axis",smartA$axis1)
      if (percent & is.null(smartA$axis1$labels))
        smartA$axis2$labels <- paste(100*smartA$axis2$at,"%")
      do.call("axis",smartA$axis2)
    }
  }

  # }}}
  # {{{ adding lines

  nix <- lapply(1:nlines, function(s) {
    lines(x = X, y = y[s,], type = "s", col = col[s], lty = lty[s], lwd = lwd[s])
  })

  # }}}
  # {{{ legend

  if(legend==TRUE && !add && !is.null(rownames(y))){
    save.xpd <- par()$xpd
    do.call("legend",smartA$legend)
    par(xpd=save.xpd)
  }
  
  # }}}

  invisible(x)
}

