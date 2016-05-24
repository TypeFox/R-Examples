# does a nice plot of two density curves with overlap shaded

overlapPlot <-
function(A, B, xscale=24, xcenter=c("noon", "midnight"),
    linetype=c(1, 2), linecol=c('black', 'blue'),
    linewidth=c(1,1), olapcol='lightgrey', rug=FALSE, extend=NULL,
    n.grid=128, kmax = 3, adjust = 1, ...)  {
    # xlab="Time", ylab="Density", ylim, now passed via "..."

  isMidnt <- match.arg(xcenter) == "midnight"

  bwA <- getBandWidth(A, kmax=kmax) / adjust
  bwB <- getBandWidth(B, kmax=kmax) / adjust
  if(is.na(bwA) || is.na(bwB))
    stop("Bandwidth estimation failed.")
  xsc <- if(is.na(xscale)) 1 else xscale / (2*pi)
  # xxRad <- seq(0, 2*pi, length=n.grid)
  if (is.null(extend)) {
    xxRad <- seq(0, 2*pi, length=n.grid)
  } else {
    xxRad <- seq(-pi/4, 9*pi/4, length=n.grid)
  }
  if(isMidnt)
    xxRad <- xxRad - pi
  xx <- xxRad * xsc
  densA <- densityFit(A, xxRad, bwA) / xsc
  densB <- densityFit(B, xxRad, bwB) / xsc
  densOL <- pmin(densA, densB)

  # Deal with ... argument:
  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]
  defaultArgs <- list(
    main=paste(deparse(substitute(A)), "and", deparse(substitute(B))), 
    xlab="Time", ylab="Density",
    bty="o", type="l", xlim=range(xx), ylim = c(0, max(densA, densB)))
  useArgs <- modifyList(defaultArgs, dots)

  selPlot <- names(useArgs) %in%
    c(names(as.list(args(plot.default))), names(par(no.readonly=TRUE)))
  plotArgs <- useArgs[selPlot]
  plotArgs$x <- 0
  plotArgs$y <- 0
  plotArgs$type <- "n"
  plotArgs$xaxt <- "n"
  do.call(plot, plotArgs, quote=TRUE)

  plotTimeAxis(xscale, ...)
  polygon(c(max(xx), min(xx), xx), c(0, 0, densOL), border=NA, col=olapcol)
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

  # if(rug)
  segments(xx[1], 0, xx[n.grid], 0, lwd=0.5)
  lines(xx, densA, lty=linetype[1], col=linecol[1], lwd=linewidth[1])
  lines(xx, densB, lty=linetype[2], col=linecol[2], lwd=linewidth[2])
  if(rug) {
    if(isMidnt) {
      A <- ifelse(A < pi, A, A - 2*pi)
      B <- ifelse(B < pi, B, B - 2*pi)
    }
    axis(1, at=A*xsc, labels=FALSE, tcl= 0.35, lwd=0, lwd.ticks=0.5, col=linecol[1])
    axis(1, at=B*xsc, labels=FALSE, tcl=-0.35, lwd=0, lwd.ticks=0.5, pos=0, col=linecol[2])
  }
  return(invisible(data.frame(x = xx, densityA = densA, densityB = densB)))
}
