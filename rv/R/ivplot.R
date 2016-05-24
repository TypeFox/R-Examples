.plot_interval <- function (data, Intervals, use_color=TRUE, xlim, lwdex=1, cexex=1)
{
  # Intervals
  data$.use_color <- use_color
  data$.lwdex <- lwdex
  data$.cexex <- cexex
  for (A in Intervals) {
    x_columns <- A$x
    D <- do.call(.uniquelist, c(A, data))
    x <- unlist(data[x_columns], use.names=FALSE)
    if (any(is.na(x))) next
    D$x <- .fix_x(x, xlim=xlim)
    args <- do.call(.interval_args, D)
    do.call(do.call, args=args)
  }
}

.prepDevice <- function (name, device="X11", hwx=1.5, size=c(210,297)/25.4/2, font.family="Courier")
{
  hw <- (size*hwx)
  if (device=="pdf") {
    if (regexpr("\\.pdf$", name) > 0) {
      filename <- name
    } else {
      filename <- paste(name, ".pdf", sep="")
    }
    pdf(
      file   = filename,
      height = hw[1],
      width  = hw[2],
      onefile = TRUE, 
      family = font.family,
      title  = "",
     paper  = "special"
    )
  } else if (device=="X11") {
    filename <- ""
    X11(
      height = hw[1],
      width  = hw[2],
      pointsize = 12
    )
  } else {
    stop("no such device:", device)
  }
  return(filename)
}

.uniquelist <- function (...)
{
  A <- list(...)
  L <- list()
  for (name in unique(names(A))) {
    a <- A[[name]]
    names(a) <- NULL
    L[[name]] <- A[[name]]
  }
  return(L)
}

.interval_args <- function (x, y, color="black", lty="solid", symbol=19, lwd=1, cex=1, .use_color=TRUE, bg=color, .lwdex=1, .cexex=1, ...)
{
  A <- list(col = if (.use_color) color else "black")
  a_point <- (length(x)==1)
  A <- if (a_point) {
    list(what="points", args=c(A, pch=symbol, x=x[1], y=y[1], cex=cex*.cexex, bg=bg))
  } else {
    list(what="segments", args=c(A, lty=lty, lwd=lwd*.lwdex, x0=x[1], x1=x[2], y0=y[1], y1=y[1]))
  }
  return(A)
}

.fix_x <- function (x, xlim, f=0.0125)
{
  if (length(x)==1) return(x)
  x <- sort(x)
#print(x);browser()
  d <- (xlim[2]-xlim[1])*f
  x[is.infinite(x)] <- xlim[is.infinite(x)]
  if (x[2]<xlim[1]) x[2] <- (xlim[1]+d)
  if (x[1]>xlim[2]) x[1] <- (xlim[2]-d)
  return(x)
}


ivplot <- function (X, name="", file.name="", split=NULL, Intervals=NULL, xlim, left.margin=3, x.ticks=NULL, exp.labels=FALSE, xlab="", title="", top.axis=FALSE, use_color=TRUE, vline=NULL, device="X11", size=c(297,210)/25.4/2, font.family="Courier", cex.label=NULL, ...)
{
# X : data frame or a list of data frames
# split : name of the column by which to split the data frame
# Intervals : list of lists specifying intervals to output e.g. list(x=c("Q2.5", "Q97.5")).
#
  .hilo <- function (high=1, low=0.5, n=nrow(X)) {
    min(high, max(low, high-((n-20)/80)*(high-low)))
  }
  if (is.list(X) && all(sapply(X, is.data.frame))) {
    A <- X
  } else if (is.data.frame(X)) {
    if (is.null(split)) {
      split <- names(X)[1]
    }
    if (!(is.character(split) && all(split%in%names(X)))) {
      stop()
    }
    if (length(X)==0) {
      stop()
    }
    A <- split(X, X[[split[1]]])
  } else {
    stop("X must be a data frame or a list of data frames")
  }
  X <- do.call(rbind, A) ## Re-join to ensure right order!
  n_levels      <- length(A)
  y_labels      <- names(A)
  y_lengths     <- sapply(A, nrow)
  y_end_coord   <- cumsum(y_lengths)
  y_start_coord <- (1 + c(0, y_end_coord[-n_levels]))
  #
  #
  gap <- (2.5-.hilo(2, 0))
  if (is.null(X$y)) {
    X$y <- seq(from=1, to=nrow(X))
    shift <- unlist(lapply(seq_along(y_lengths), function (i) rep((i-1)*gap, y_lengths[i])))
    X$y <- (X$y + shift)
    y_start_coord <- X$y[y_start_coord]
  }
  ylim <- range(X$y)
  #
  ylim[2] <- (ylim[2]+1)
  #
  file.name <- .prepDevice(name=file.name, device=device, size=size, font.family=font.family)
  #
  old_mar <- par("mar")
  new_mar <- c(5.1, left.margin+5.1, 2.1, 1.1)
  #
  old_par <- par(mar = new_mar)
  #
  .cleanup <- function () {
    par(old_par)
    if (device=="pdf") dev.off()
  }
  on.exit(.cleanup())
  #
  las <- 1
  #
  all_columns <- unlist(lapply(Intervals, "[[", "x"))
  if (!all(all_columns%in%names(X))) {
    stop("The argument `Intervals' is invalid")
  }
  all_x <- unlist(X[all_columns])
  all_finite_x <- all_x[is.finite(all_x)]
  range_all_finite_x <- range(all_finite_x)
  if (missing(xlim)) { 
    xlim <- range_all_finite_x
  } else {
    xlim <- xlim[1:2]
    if (any(is.na(xlim))) {
      xlim[is.na(xlim)] <- range_all_finite_x[is.na(xlim)]
    }
  }
  x_row_coords <- if (is.null(x.ticks)) pretty(xlim) else x.ticks
  x_label_num <- signif(if (exp.labels) exp(x_row_coords) else x_row_coords, 3)
  if (!is.null(x.ticks)) {
    x_labels <- if (length(names(x.ticks))>0) names(x.ticks) else paste(x_label_num)
  } else {
    x_labels <- paste(x_label_num)
  }
  #
  Dots <- list(...)
  Args <- .uniquelist(..., 
    x=xlim, y=ylim, xlim=xlim, ylim=rev(ylim), 
    las=las, type="n", axes=FALSE, xlab=xlab, ylab=""
  )
  if ("cex.main" %in% names(Args)) {
    cex.main <- Args$cex.main
    Args$cex.main <- NULL
    Dots$cex.main <- NULL
  } else {
    cex.main <- (1.0)
  }
  if ("cex.lab" %in% names(Args)) {
    Dots$cex.lab <- NULL
  }
  do.call(plot, Args)
  #
  # box()
  axis(1, at = x_row_coords, labels=x_labels, cex.axis=0.66)
  if (top.axis) {
    axis(3, at = x_row_coords, labels=x_labels, cex.axis=0.66)
  }
  #axis(2, at = y_start_coord, labels = y_labels, tick = FALSE, 
  #  line = FALSE, pos = NA, outer = FALSE, font = NA, 
  #  hadj = 0, las = 1)
  #
  for (i in seq_along(X)) {
    if (is.factor(X[[i]])) {
      X[[i]] <- as.character(X[[i]])
    }
  }
  cex1 <- .hilo(1, 0.8)
  cex2 <- .hilo(1, 0.45)
  if (! is.null(cex.label)) {
    cex2 <- cex2 * cex.label
  }
  lwdex <- .hilo(1.25, 0.5)
  cexex <- .hilo(1.25, 0.5)
  labels1 <- y_labels
  labels1.y <- (y_start_coord - 0.5)
  title(main=title, cex.main=cex.main)
  if (length(X$label)>0) {
    mtext(text=X$label, side=2, line=0, at=X$y, las=1, adj=1, cex=cex2, font=1)
  }
  mtext(text=labels1, side=2, line=left.margin + 0L, at=labels1.y, las=1, padj=0, adj=1, cex=cex1, font=2)
  if (length(vline)>0) {
    for (i in seq_along(vline)) {
      v <- vline[i]
      lty <- names(vline)[i]
      if (is.null(lty) || is.na(lty)) {
        lty <- "dotted"
      }
      abline(v=v, lty=lty, lwd=0.50)
    }
  }
  A <- Dots
  A$Intervals <- Intervals
  A$use_color <- use_color
  A$xlim <- xlim
  A$lwdex <- lwdex
  A$cexex <- cexex
  for (i in 1:nrow(X)) {
    A$data <- as.list(X[i,])
    do.call(.plot_interval, A)
  }
  if (length(X$label2)>0) {
    mtext(text=X$label2, side=2, line=-0.25, at=X$y, las=1, adj=0, cex=cex2 * 0.9, font=1)
  }
  invisible(file.name)
}
