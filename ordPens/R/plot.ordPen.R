plot.ordPen <- function(x, whichlam = NULL, whichx = NULL, type = NULL,
xlab = NULL, ylab = NULL, main = NULL, xlim = NULL, ylim = NULL, col = NULL, ...)
  {
    px <- length(x$xlevels)
    xgrp <- rep(1:px,x$xlevels)
    tol <- .Machine$double.eps^0.5
    if (is.null(whichlam))
      whichlam <- 1:ncol(x$coef)
    else if (!is.numeric(whichlam) | max(whichlam) > ncol(x$coef) |
    any(abs(whichlam - round(whichlam)) > tol))
      stop("incorrect whichlam")

    if (rownames(x$coef)[1] == "intercept")
      xcoefs <- x$coef[2:(length(xgrp)+1),whichlam,drop=FALSE]
    else
      xcoefs <- x$coef[1:length(xgrp),whichlam,drop=FALSE]

    if (is.null(whichx))
      whichx <- 1:px
    else if (!is.numeric(whichx) | max(whichx) > px |
    any(abs(whichx - round(whichx)) > tol))
      stop("incorrect whichx")

    if (is.null(xlab))
      xlab <- "level"

    if (is.null(ylab))
      ylab <- "dummy coefficient"

    if (is.null(type))
      type <- "b"
     
    noylims <- is.null(ylim)
    nomain <- is.null(main)
    nocol <- is.null(col)
    multcol <- length(col) > 1
      
    if (nocol)
      cols <- grey(seq(0,0.7,length=length(whichlam)))
    else if (multcol)
      {
        if (length(col) != length(whichlam))
          stop("incorrect length(col)")
        else
          cols <- col
      }
      
    devAskNewPage(length(whichx)>1)
    for (wx in whichx)
      {
        xlam <- xcoefs[xgrp==wx, ,drop=FALSE]
        
        if (noylims)
          ylim <- c(min(xlam),max(xlam))
          
        if (nomain)
          {
            xname <- rownames(xlam)[1]
            main <- substr(xname,1,nchar(xname)-2)
          }

        if (nocol | multcol)
          col <- cols[1]

        plot(1:nrow(xlam), xlam[,1], xlim = xlim, ylim = ylim, main = main,
        xlab = xlab, ylab = ylab, type = type, col = col, ...)
        if (ncol(xlam) > 1)
          {
            for (wl in 2:ncol(xlam))
              {
                if (nocol | multcol)
                  col <- cols[wl]
                lines(1:nrow(xlam), xlam[,wl], type = type, col = col, ...)
              }
          }
      }
  }