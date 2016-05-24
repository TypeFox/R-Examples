# Plots a kernel density for circular data

# A: a sample of times of observations in radians
# adjust: smoothing parameter (adjust = 1/c in old code)

densityPlot <-
function(A, xscale=24, xcenter=c("noon", "midnight"),
    add=FALSE, rug=FALSE, extend="lightgrey",
    n.grid=128, kmax = 3, adjust = 1, ...)  {
    # ylim, xlab="Time", ylab="Density" now included in defaultArgs

  isMidnt <- match.arg(xcenter) == "midnight"
  bw <- getBandWidth(A, kmax=kmax) / adjust
  if(is.na(bw))
    stop("Bandwidth estimation failed.")
  # xx <- seq(0, 2*pi, length=n.grid)
  if (is.null(extend)) {
    xx <- seq(0, 2*pi, length=n.grid)
  } else {
    xx <- seq(-pi/4, 9*pi/4, length=n.grid)
  }
  if(isMidnt)
    xx <- xx - pi
  densA <- densityFit(A, xx, bw)
  xsc <- if(is.na(xscale)) 1 else xscale / (2*pi)
  toPlot <- cbind(x = xx * xsc, y = densA / xsc)

  # Deal with ... argument:
  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]
  defaultArgs <- list(main=deparse(substitute(A)), bty="o", type="l",
    xlab="Time", ylab="Density", ylim = c(0, max(toPlot[,'y'])))
  useArgs <- modifyList(defaultArgs, dots)

  if(!add)  {
    selPlot <- names(useArgs) %in%
      c(names(as.list(args(plot.default))), names(par(no.readonly=TRUE)))
    plotArgs <- useArgs[selPlot]
    plotArgs$x <- toPlot
    plotArgs$y <- NULL
    plotArgs$type <- "n"
    plotArgs$xaxt <- "n"
    do.call(plot, plotArgs, quote=TRUE)

    plotTimeAxis(xscale, ...)
    abline(h=0, col='grey')
    if(!is.null(extend)) {
      if(isMidnt) {
        wrap <- c(-pi, pi) * xsc
      } else {
        wrap <- c(0, 2*pi) * xsc
      }
      edge <- par('usr')
      rect(c(edge[1], wrap[2]), rep(edge[3], 2), c(wrap[1], edge[2]), rep(edge[4],2),
        border=NA, col=extend)
      box(bty=useArgs$bty)
    }
  }
  selPlot <- names(useArgs) %in% names(par(no.readonly=TRUE))
  plotArgs <- useArgs[selPlot]
  plotArgs$x <- toPlot
  plotArgs$y <- NULL
  do.call(lines, plotArgs, quote=TRUE)

  if(rug)  {
    if(isMidnt)
      A <- ifelse(A < pi, A, A - 2*pi)
    rug(A * xsc, ...)  # do.call(rug, doesn't work !!
  }
  return(invisible(as.data.frame(toPlot)))
}
