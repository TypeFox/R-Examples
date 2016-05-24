###########################################################################
##                                                                       ##
## screeplot     - 'screeplot' methods.                                  ##
##                                                                       ##
## Created       : 27-May-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1                                                   ##
## Last modified : 27-May-2006                                           ##
##                                                                       ##
## ARGUMENTS:                                                            ##
## x             - object on which method dispatch applied               ##
## k             - number of analogues to use. If missing 'k' is chosen  ##
##                 automatically as the 'k' that achieves lowest RMSE.   ##
## restrict      - logical; restrict comparison of k-closest model to    ##
##                 k <= 'restrict'.                                      ##
## display       - which aspect of 'mat' object to plot. One of 'rmse',  ##
##                 'rmsep', 'avg.bias', 'max.bias', or 'r.squared'.      ##
## weighted      - Logical. Should the analysis use weighted mean of env ##
##                 data of analogues as fitted/estimated values?         ##
## xlab, ylab    - x- and y-axis labels                                  ##
## main, sub     - main and subtitle for the plot                        ##
## ...           - arguments passed to other graphics functions          ##
##                                                                       ##
###########################################################################
screeplot.mat <- function(x,
                          k,
                          restrict = 20,
                          display = c("rmsep","avg.bias",
                            "max.bias","r.squared"),
                          weighted = FALSE,
                          col = "red",
                          xlab = NULL,
                          ylab = NULL,
                          main = NULL,
                          sub = NULL,
                          ...)
  {
    if (!inherits(x, "mat")) 
      stop("use only with \"mat\" objects")
    if(missing(k) || is.null(k)) {
      n.obs <- nrow(x$orig.x)
      if (n.obs > restrict)
        k <- restrict
      else
        k <- n.obs
    }
    if(missing(display))
      display <- "rmsep"
    display <- match.arg(display)
    captions <- c("RMSEP","Average bias","Maximum bias","R squared")
    names(captions) <- c("rmsep","avg.bias","max.bias","r.squared")
    if (is.null(xlab))
      xlab <- "No. of analogues, k"
    if (is.null(ylab))
      ylab <- captions[display]
    if (is.null(main))
      main <- deparse(substitute(x))
    if(is.null(sub)) {
      cal <- x$call
      if (!is.na(m.f <- match("formula", names(cal)))) {
        cal <- cal[c(1, m.f)]
        names(cal)[2] <- ""
      }
      cc <- deparse(cal, 80)
      nc <- nchar(cc[1])
      abbr <- length(cc) > 1 || nc > 75
      sub <- if (abbr) 
        paste(substr(cc[1], 1, min(75, nc)), "...")
      else cc[1]
    }
    #if(display == "rmsep")
    #  display <- "rmsep"
    if(weighted)
      dat <- x$weighted[[display]][1:k]
    else
      dat <- x$standard[[display]][1:k]
    plot(1:k, dat, type = "n", ylab = ylab, xlab = xlab, main = main,
         sub = sub, ...)
    if(restrict > 20) {
      lines(1:k, dat, type = "l", ...)
    } else {
      lines(1:k, dat, type = "b", pch = "", ...)
      text(1:k, dat, labels = as.character(seq(1, k)), cex = 0.8, ...)
    }
    invisible()
  }

##screeplot.bootstrap <- function(x,
screeplot.bootstrap.mat <- function(x,
                                    k,
                                    restrict = 20,
                                    display = c("rmsep",
                                      "avg.bias","max.bias","r.squared"),
                                    legend = TRUE,
                                    loc.legend = "topright",
                                    col = c("red", "blue"),
                                    xlab = NULL,
                                    ylab = NULL,
                                    main = NULL,
                                    sub = NULL,
                                    ..., 
                                    lty = c("solid","dashed")
                                    ) {
  if (!inherits(x, "bootstrap.mat")) 
    stop("use only with \"bootstrap.mat\" objects")
  ##if (!inherits(x, "bootstrap")) 
  ##  stop("use only with \"bootstrap\" objects")
  if(missing(k) || is.null(k)) {
    n.obs <- length(x$model$rmsep)
    if (n.obs > restrict)
      k <- restrict
    else
      k <- n.obs
  }
  ##dotargs <- list(...)
  if(missing(display))
    display <- "rmsep"
  else
    display <- match.arg(display)
  captions <- c("Error","Average bias","Maximum bias","R squared")
  names(captions) <- c("rmsep","avg.bias","max.bias","r.squared")
  if(length(col) == 1)
    col <- rep(col, 2)
  if (is.null(xlab))
    xlab <- "No. of analogues, k"
  if (is.null(ylab))
    ylab <- captions[display]
  if (is.null(main))
    main <- deparse(substitute(x))
  if(is.null(sub)) {
    cal <- x$call
    if (!is.na(m.f <- match("formula", names(cal)))) {
      cal <- cal[c(1, m.f)]
      names(cal)[2] <- ""
    }
    cc <- deparse(cal, 80)
    nc <- nchar(cc[1])
    abbr <- length(cc) > 1 || nc > 75
    sub <- if (abbr) 
      paste(substr(cc[1], 1, min(75, nc)), "...")
    else cc[1]
  }
  if(display == "rmsep") {
    dat1 <- x$model$rmsep[1:k]
    dat2 <- x$bootstrap$rmsep[1:k]
  } else {
    dat1 <- x$model[[display]][1:k]
    dat2 <- x$bootstrap[[display]][1:k]
  }
  ylims <- range(dat1, dat2)
  ## add 10% of range if legend
  if(legend)
    ylims[2] <- ylims[2] + ((ylims[2] - ylims[1]) * 0.1)
  plot(1:k, dat1, type = "n", ylim = ylims, ylab = ylab, xlab = xlab,
       main = main, sub = sub, col = col[1], ...)
  if(restrict > 20) {
    lines(1:k, dat1, type = "l", lty = lty[1], col = col[1], ...)
    lines(1:k, dat2, type = "l", lty = lty[2], col = col[1], ...)
  } else {
    lines(1:k, dat1, type = "b", pch = "", lty = lty[1], col = col[1], ...)
    text(1:k, dat1, labels = as.character(seq(1, k)), cex = 0.8, ...)
    lines(1:k, dat2, type = "b", pch = "", lty = lty[2], col = col[2], ...)
    text(1:k, dat2, labels = as.character(seq(1, k)), cex = 0.8, ...)
  }
  if(legend) {
    legend(loc.legend, legend = c("LOO","Bootstrap"), lty = lty,
           inset = 0.01, bty = "n", horiz = FALSE, cex = 0.8,
           col = col)
  }
  invisible()
}
