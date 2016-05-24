"plot.seas.sum" <-
function(x, var, norm = "days", year.filter, ylim,
         start = 1, rep = 0, col = "lightgrey", main, ylab, ...) {
  orig <- as.character(substitute(x))[[1]]
  x <- seas.sum.check(x, orig, var, norm, year.filter)
  var <- x$var
  val.m <- x$seas[,,var]/x$norm[,,var]  # normalize and deconstruct array
  val.m[!is.finite(val.m)] <- NA
  val.d <- data.frame(val.m)
  names(val.d) <- x$bins
  val <- stack(val.d)
  val$ind <- factor(val$ind, levels=x$bins)  # correct order
  num <- length(x$bins)
  if (missing(main)) {
    fun <- sprintf("%s/%s", var, norm)
    main <- .seastitle(id=x$id, name=x$name, orig=orig, fun=fun,
                       range=x$year.range)
  }
  xlab <- .seasxlab(x$width, x$start.day)
  units <- if (is.null(x$units[[var]])) NULL else gettextf("%s/day", x$units[[var]])
  if (missing(ylab))
    ylab <- .seasylab(var, long.name=x$long.name[[var]],units)
  #if (add.alt) {
  #  mar <- c(5.1, 4.1, 4.1, 4.1)
  #  bty <- "u"
  #  alt.ylab <- gettextf("%s/day", alt.unit)
  #} else {
  #mar <- c(5.1, 4.1, 4.1, 2.1)
  #bty <- "l"
  #}
  log <- if (par("ylog")) "y" else ""
  if (missing(ylim))
    ylim <- range(c(0, val$values), na.rm=TRUE) * 1.04
  else if (length(ylim) == 1)
    ylim <- c(0, ylim)
  if (log=="y")
    ylim[1] <- min(val$values, na.rm=TRUE)
  par(xaxs="r", yaxt="n", xaxt="n")
  frame()
  plot.window(xlim=c(0.5, num + 0.5), ylim=ylim, log=log)
  if (start != 1) {
    warning("'start' not yet supported; must be 1 for now")
    start <- 1
  }
  if (rep != 0) {
    warning("'rep' not yet supported; must be 0 for now")
    rep <- 0
  }
  .seasmonthgrid(x$width, x$bin.lengths, start, rep, x$start.day)
  op <- getOption("seas.bxp")
  pl <- boxplot(values ~ ind, val, col=col, varwidth=TRUE, names=NA,
                add=TRUE, outcex=op$outcex * par("cex"))
  #if (add.alt) {
  #  alt.ax <- pretty(ylim * slope + inter)
  #  axis(side=4, at=(alt.ax - inter) / slope, lab=alt.ax, srt=90)
  #  if (!is.na(ylab[2]))
  #    mtext(ylab[2], side=4, line=2.8)
  #}
  par(yaxt="s", xaxt="s")
  axis(2, lwd=par("lwd"))
  axis(1, 1:num, x$bins, lwd=par("lwd"))
  title(main, xlab=xlab, ylab=ylab[1])
  invisible(pl)
}
