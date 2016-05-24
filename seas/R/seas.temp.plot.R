"seas.temp.plot" <-
function(x, width = 11, start = 1, rep = 0, start.day = 1,
         var = c("t_min", "t_max", "t_mean"),
         add.alt = FALSE, ylim, main, ylab, ...) {
  orig <- as.character(substitute(x))[[1]]
  sc <- seas.df.check(x, orig, var[1:2])  # only need to check first two vars
  op <- getOption("seas.temp")
  if (length(var) == 1) {
    stop("need at least 't_min' and 't_max' variables")
  } else if (length(var) == 2) {
    warning("calculating mean temperatures from min and max temperatures")
    x$t_mean <- rowMeans(x[,var[1:2]])
    var[3] <- "t_mean"
  } else {
    x$t_mean <- x[[var[3]]]
  }
  seas.df.check(x, orig, var[3])
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
  if (missing(main))
    main <- sc$main
  xlab <- .seasxlab(width, start.day)
  at.label <- ((1:num) + start - 2) %% num.fact + 1
  if (!is.numeric(width))
    at.label <- levels(x$fact)[at.label]
  degSymb <- iconv("\260", "latin1", "")
  if (is.null(sc$units))  # assume degC
    sc$units <- sprintf("%sC", degSymb)
  if (is.null(sc$long.name))
    sc$long.name <- "temperature"
  unit <- substring(sc$units, nchar(sc$units, type="chars"))
  if (!unit %in% c("C", "F", "K"))
    warning(
      paste(
        gettextf("the unit %s in %s is not recognized",
                 sQuote(sc$units),
                 sQuote(sprintf("attr(%s$%s, \"units\")",
                 orig, var[1]))),
        gettextf("this must end with %s, %s or %s",
                 sQuote("C"), sQuote("F"), sQuote("K")), sep="\n ... "))
  if (missing(ylab))
    ylab <- .seasylab(orig, "Temperature", sc$units)
  mar <- par("mar")
  ylog <- par("ylog")
  if (add.alt) {
    mar[4] <- mar[2]
    bty <- "u"
  } else {
    bty <- "l"
  }
  par(mar=mar, bty=bty, xaxs="r", yaxt="n", xaxt="n")
  if (missing(ylim))
    ylim <- range(x$t_mean, na.rm=TRUE)
  frame()
  plot.window(xlim=c(0.5, num + 0.5), ylim)
  .seasmonthgrid(width, days, start, rep, start.day)
  varwidth <- TRUE
  lwd <- par("lwd")
  seas.bxp <- function(at) {
    pl <- boxplot(t_mean ~ fact, x, at=at, col=op$col[1],
                  log=ylog, varwidth=varwidth, names=NA, add=TRUE,
                  outcex=getOption("seas.bxp")$outcex * par("cex"))
    # compute mean diurnal variability
    dmin <- tapply(x[,var[1]], x$fact, mean, na.rm=TRUE)
    dmax <- tapply(x[,var[2]], x$fact, mean, na.rm=TRUE)
    segments(at, dmax, at, dmin, col=op$col[2],
             lwd=op$lwd * lwd, lend="square")
    invisible(pl)
  }
  pl <- seas.bxp(1:num.fact)
  if (rep>0) {
    n <- 1
    if (rep >= num.fact) {  # whole repetitions
      for (n in 1:floor(rep / num.fact))
        seas.bxp((num.fact * n + 1):(num.fact * (n + 1)))
      n <- n + 1
    }
    if (rep %% num.fact != 0) {  # non-whole repetitions
      lv <- levels(x$fact)[1:(rep %% num.fact)]
      x <- x[x$fact %in% lv,]
      x$fact <- factor(x$fact, lv)
      seas.bxp((num.fact * n + 1):(num.fact * n + rep %% num.fact))
    }
  }
  par(yaxt="s", xaxt="s")
  axis(2, lwd=lwd)
  axis(1, 1:num, at.label, lwd=lwd)
  title(main, xlab=xlab, ylab=ylab)
  if (unit == "C")
    abline(h=0)
  else if (unit == "F")
    abline(h=32)
  else if (unit == "K")
    abline(h=273.15)
  if (add.alt) {
    f2c <- function(v)((v - 32) / 1.8)
    c2f <- function(v)(1.8 * v + 32)
    if (unit == "C") {
      alt.ax <- pretty(c2f(ylim))
      alt.at <- f2c(alt.ax)
      alt.unit <- sprintf("%sF", degSymb)
    } else if (unit == "F") {
      alt.ax <- pretty(f2c(ylim))
      alt.at <- c2f(alt.ax)
      alt.unit <- sprintf("%sC", degSymb)
    } else {  # degK
      alt.ax <- pretty(ylim - 273.15)
      alt.at <- alt.ax + 273.15
      alt.unit <- sprintf("%sC", degSymb)
    }
    axis(side=4, at=alt.at, labels=alt.ax, srt=90, lwd=lwd)
    mtext(alt.unit, side=4, line=par("mgp")[1])
  }
  invisible(pl)
}
