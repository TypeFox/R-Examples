"seas.var.plot" <-
function(x, var, width = 11, start = 1, rep = 0, start.day = 1,
         col, ylim, add.alt, alt.ylab, main, ylab, ylog, ...) {
  orig <- as.character(substitute(x))[[1]]
  if (missing(var))
    stop(gettextf("%s must be supplied", sQuote("var")))
  sc <- seas.df.check(x, orig, var)
  x$var <- x[[var[1]]]
  if (missing(main))
    main <- sc$main
  x$fact <- mkseas(x, width, start.day)
  days <- attr(x$fact, "bin.lengths")
  num.fact <- length(levels(x$fact))
  num <- num.fact + rep
  if (start < 1 || start >= num.fact) {
    warning("'start' should be >= 1 and < to the number of bins")
    start <- 1
  }
  if (start > 1)  # re-order
    x$fact <- factor(x$fact, levels(x$fact)[(1:num.fact - 2 + start) %% num.fact + 1])
  xlab <- .seasxlab(width, start.day)
  at.label <- ((1:num) + start - 2) %% num.fact + 1
  if (!is.numeric(width))
    at.label <- levels(x$fact)[at.label]
  if (missing(ylog) || !is.logical(ylog))
    ylog <- FALSE
  if (ylog) {  # transform the data; different boxplot than par(ylog=TRUE)
    if (min(x$var, na.rm=TRUE) <= 0)
      warning("values of 'x' will be removed for log10 transform")
    x <- x[x$var > 0,]
    x$ovar <- x$var
    x$var <- log10(x$ovar)
  }
  if (missing(ylab))
    ylab <- sc$ylab
  if (missing(add.alt))
    add.alt <- FALSE
  else {
    if (is.numeric(add.alt) && length(add.alt) == 2) {
      slope <- add.alt[1]
      inter <- add.alt[2]
      add.alt <- TRUE
    } else if (!ylog) {
      warning(gettextf(paste("'add.alt' must give the slope and intercept",
                             "between the primary and secondary axis;",
                             " use %s\n", sep="\n"),
                       sQuote("add.alt=c(slope, inter)")))
      add.alt <- FALSE
    }
  }
  orig.par <- par(no.readonly=TRUE)
  mar <- orig.par$mar
  log <- if (orig.par$ylog) "y" else ""
  if (add.alt) {
    if (missing(alt.ylab))
      alt.ylab <- NULL
    mar[4] <- mar[2]  # copy left-hand side margin
    bty <- "u"
  } else {
    bty <- "l"
  }
  par(mar=mar, bty=bty, xaxs="r", yaxt="n", xaxt="n")
  if (missing(ylim))
    ylim <- range(x$var, na.rm=TRUE)
  else if (ylog)
    ylim <- log10(ylim)
  frame()
  if (ylog) {  # need to determine yaxp from original values
    plot.window(xlim=c(0, 1), 10^ylim, log="y")
    yaxp <- par("yaxp")
    ylax <- axTicks(2)
    par(ylog=FALSE)  # set back
  }
  plot.window(xlim=c(0.5, num + 0.5), ylim, log)
  .seasmonthgrid(width, days, start, rep, start.day)
  op <- getOption("seas.bxp")
  if (missing(col))
    col <- op$boxcol
  varwidth = TRUE
  seas.boxplot <- function(at) {
    pl <- boxplot(var ~ fact, x, at=at, col=col, varwidth=varwidth,
                  names=NA, add=TRUE, outcex=op$outcex * par("cex"))
    invisible(pl)
  }
  pl <- seas.boxplot(1:num.fact)
  if (rep > 0) {
    n <- 1
    if (rep>=num.fact) {  # whole repetitions
      for (n in 1:floor(rep / num.fact))
        seas.boxplot((num.fact * n + 1):(num.fact * (n + 1)))
      n <- n + 1
    }
    if (rep %% num.fact != 0) {  # non-whole repetitions
      lv <- levels(x$fact)[1:(rep %% num.fact)]
      x <- x[x$fact %in% lv,]
      x$fact <- factor(x$fact, lv)
      seas.boxplot((num.fact * n + 1):(num.fact * n + rep %% num.fact))
    }
  }
  par(yaxt=orig.par$yaxt, xaxt=orig.par$xaxt)
  if (ylog)
    axis(2, log10(ylax), ylax, lwd=par("lwd"))
  else
    axis(2, lwd=par("lwd"))
  axis(1, 1:num, at.label, lwd=par("lwd"))
  title(main, xlab=xlab, ylab=ylab)
  if (add.alt) {
    if (ylog) {  # done differently for ylog
      axis(4)
    } else {
      alt.ax <- pretty(ylim * slope + inter)
      axis(side=4, at=(alt.ax - inter) / slope, labels=alt.ax,
           srt=90, lwd=par("lwd"))
    }
    if (!is.null(alt.ylab))
      mtext(alt.ylab, side=4, line=par("mgp")[1])
  }
  invisible(pl)
}
