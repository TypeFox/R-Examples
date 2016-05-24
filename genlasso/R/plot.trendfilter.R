plot.trendfilter <- function(x, style=c("trend", "path"), lambda,
                             nlam, df, xlab, ylab, ...) {
  style = style[[1]]
  if (!(style %in% c("trend", "path"))) {
    stop("Invalid style, must be \"trend\" or \"path\".")
  }
  
  if (style=="path") {
    plot.genlasso(x,xlab=xlab,ylab=ylab,...)
  }
  else {
    if (!any(class(x)=="trendfilter")) {
      stop("Passed object must be of class \"trendfilter\".")
    }
    if (missing(xlab)) xlab = "Position"
    if (missing(ylab)) ylab = "Trend filtering estimate"

    # If nothing has been specified, default is to pick out
    # 10 lambda values
    if (missing(lambda) && missing(nlam) && missing(df)) nlam = 10
    co = coef.genlasso(x,lambda,nlam,df)
    if (length(co$lambda)==0) stop("Nothing to plot.")

    xvals = 1:nrow(co$beta)
    if (!is.null(x$pos)) xvals = x$pos
    else xvals = 1:nrow(co$beta)
    
    # If there's no X matrix, draw y
    if (is.null(x$X)) {
      plot(xvals, x$y, xlab=xlab, ylab=ylab, ...)
      matplot(xvals, co$beta, type="l", lty=1, add=TRUE)
    }
    # Otherwise don't draw y
    else {
      matplot(xvals, co$beta, xlab=xlab, ylab=ylab, type="l", lty=1, ...)
    }

    invisible(co)
  }
}

