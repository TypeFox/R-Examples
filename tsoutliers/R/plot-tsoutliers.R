
#see arguments x.axis and args.plot(xaxt="n"), decide which one should prevail
#
#argument ellipsis '...' is not used
#this approach (args.plot,...) is used because if plot arguments can affect 
#different parts of the plot, e.g. col="blue" for the original series or the outlier effects,...
#there may be a given argumetn e.g. "col" with different values for each element in the plot

#   args.lines.y = list(col = "gray80")
#   args.lines.yadj = list(col = "blue")
#   args.lines.effects = list(type = "s", col = "red")
#   args.points = list(col = "gray80", bg = "red", pch = 21)
#   plot.points = TRUE
#   x.axis = TRUE
#   y.axis = TRUE

plot.tsoutliers <- function(x, 
  args.lines.y = list(col = "gray80"), args.lines.yadj = list(col = "blue"),
  args.lines.effects = list(type = "s", col = "red"),   
  args.points = list(col = "gray80", bg = "red", pch = 21), plot.points = TRUE, 
  args.x.axis = list(at = pretty(time(x$y)), tcl = -0.5, lwd = 0, lwd.ticks = 1),
  args.y.axis = list(at = pretty(x$y), tcl = -0.5, lwd = 0, lwd.ticks = 1),
  args.effects.axis = list(at = pretty(x$effects), tcl = -0.5, lwd = 0, lwd.ticks = 1),
  ...)
{
  # export(plot.tsoutliers) is required in NAMESPACE
  # in order to use formals(plot.tsoutliers)

  ##NOTE 
  # this approach keeps the defaults as defined above in the formals of 
  # plot.tsoutliers() rather than the defaults in plot() or lines()

  fargs.linesy <- formals(plot.tsoutliers)$args.lines.y
  efargs.linesy <- eval(fargs.linesy)
  if (!identical(args.lines.y, efargs.linesy))
  {
    args.lines.y <- c(args.lines.y, efargs.linesy)
    id <- which(duplicated(names(args.lines.y)))
    if (length(id) > 0)
      args.lines.y <- args.lines.y[-id]
  }

  fargs.linesyadj <- formals(plot.tsoutliers)$args.lines.yadj
  efargs.linesyadj <- eval(fargs.linesyadj)
  if (!identical(args.lines.yadj, efargs.linesyadj))
  {
    args.lines.yadj <- c(args.lines.yadj, efargs.linesyadj)
    id <- which(duplicated(names(args.lines.yadj)))
    if (length(id) > 0)
      args.lines.yadj <- args.lines.yadj[-id]
  }

  fargs.linesef <- formals(plot.tsoutliers)$args.lines.effects
  efargs.linesef <- eval(fargs.linesef)
  if (!identical(args.lines.effects, efargs.linesef))
  {
    args.lines.effects <- c(args.lines.effects, efargs.linesef)
    id <- which(duplicated(names(args.lines.effects)))
    if (length(id) > 0)
      args.lines.effects <- args.lines.effects[-id]
  }

  fargs.points <- formals(plot.tsoutliers)$args.points
  efargs.points <- eval(fargs.points)
  if (!identical(args.points, efargs.points))
  {
    args.points <- c(args.points, efargs.points)
    id <- which(duplicated(names(args.points)))
    if (length(id) > 0)
      args.points <- args.points[-id]
  }

  fargs.xaxis <- formals(plot.tsoutliers)$args.x.axis
  efargs.xaxis <- eval(fargs.xaxis)
  if (!identical(args.x.axis, efargs.xaxis))
  {
    args.x.axis <- c(args.x.axis, efargs.xaxis)
    id <- which(duplicated(names(args.x.axis)))
    if (length(id) > 0)
      args.x.axis <- args.x.axis[-id]
  }
  if (is.null(args.x.axis$labels))
    args.x.axis$labels <- args.x.axis$at
  args.x.axis$side <- 1

  fargs.yaxis <- formals(plot.tsoutliers)$args.y.axis
  efargs.yaxis <- eval(fargs.yaxis)
  if (!identical(args.y.axis, efargs.yaxis))
  {
    args.y.axis <- c(args.y.axis, efargs.yaxis)
    id <- which(duplicated(names(args.y.axis)))
    if (length(id) > 0)
      args.y.axis <- args.y.axis[-id]
  }
  if (is.null(args.y.axis$labels))
    args.y.axis$labels <- args.y.axis$at
  args.y.axis$side <- 2

  fargs.eaxis <- formals(plot.tsoutliers)$args.effects.axis
  efargs.eaxis <- eval(fargs.eaxis)
  if (!identical(args.effects.axis, efargs.eaxis))
  {
    args.effects.axis <- c(args.effects.axis, efargs.eaxis)
    id <- which(duplicated(names(args.effects.axis)))
    if (length(id) > 0)
      args.effects.axis <- args.effects.axis[-id]
  }
  if (is.null(args.effects.axis$labels))
    args.effects.axis$labels <- args.effects.axis$at
  if (is.null(args.effects.axis$side))
    args.effects.axis$side <- 4

  if (nrow(x$outliers) == 0)
  {
    cat(paste(sQuote("x"), "does not contain outliers to display\n"))
    return()
  }

  oldpar <- par(mar = c(0, 3, 0, 2.1), oma =  c(3, 0, 3, 0), mfcol = c(2, 1), ...)
  on.exit(par(oldpar))

  #do.call("plot", args = c(list(x = cbind(x$y, x$adj)), args.plot))
  plot(cbind(x$y, x$yadj), plot.type ="single", 
    type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  mtext(side = 3, text = "Original and adjusted series", adj = 0)

  do.call("lines", args = c(list(x = x$y), args.lines.y))
  do.call("lines", args = c(list(x = x$yadj), args.lines.yadj))

  #if (y.axis)
  #{
  #  ay <- pretty(x$y)
  #  axis(side = 2, at = ay, labels = FALSE, tcl = 0.25, lwd = 0, lwd.ticks = 1)
  #  axis(side = 2, at = ay, labels = ay, lwd = 0, lwd.ticks = 0, line = -0.5)
  #  axis(side = 2, at = ay-(ay[2] - ay[1])/2, labels = FALSE, tcl = 0.15, lwd = 0, lwd.ticks = 1)
  #}
  do.call("axis", args = args.y.axis)

  if (plot.points)
  {
    do.call("points", args = c(list(x = x$times, y = x$y[x$outliers[,"ind"]]), 
      args.points))
  }

  # bty = "u" is necessary to avoid the horizontal line between both plots 
  # looks thicker because of overlap between the lower side of top plot's box
  # and upper side of bottom plot's box

  plot(x$effects, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "u")
  do.call("lines", args = c(list(x = x$effects), args.lines.effects))
  mtext(side = 3, text = "Outlier effects", adj = 0, line = -1)

  ##NOTE
  #by default the axis are created separately with option "lwd=0" because 
  #otherwise it is noticed that there are two overlapping lines,
  #the box and the axis lines

  #if (y.axis)
  #{
  #  ay <- pretty(x$effects)
  #  axis(side = 4, at = ay, labels = FALSE, tcl = 0.25, lwd = 0, lwd.ticks = 1)
  #  axis(side = 4, at = ay, labels = ay, lwd = 0, lwd.ticks = 0, line = -0.5)
  #  axis(side = 4, at = ay-(ay[2] - ay[1])/2, labels = FALSE, tcl = 0.15, lwd = 0, lwd.ticks = 1)
  #}
  do.call("axis", args = args.effects.axis)

  #if (x.axis)
  #{
  #  ax <- pretty(time(x$y))
  #  axis(side = 1, at = ax, labels = FALSE, tcl = 0.5, lwd = 0, lwd.ticks = 1)
  #  axis(side = 1, at = ax - (ax[2] - ax[1])/2, labels = FALSE, tcl = -0.15, lwd = 0, lwd.ticks = 1)
  #  axis(side = 1, at = ax, labels = ax, lwd = 0, lwd.ticks = 0, line = -0.5)
  #}
  do.call("axis", args = args.x.axis)
}
