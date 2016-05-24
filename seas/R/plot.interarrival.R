"plot.interarrival" <-
function(x, width = 11, start = 1, rep = 0, start.day = 1,
         ylog = FALSE, maxy, main, ...) {
  orig <- as.character(substitute(x))[[1]]
  if (!inherits(x, "interarrival"))
    stop(gettextf("%s is not an %s object",
                  sQuote(orig), sQuote("interarrival")))
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  xlab <- .seasxlab(width, start.day)
  ylab1 <- gettext("wet days"); ylab2 <- gettext("dry days")
  if (missing(main))
    main <- .seastitle(id=attr(x, "id"), name=attr(x, "name"),
                       orig=orig, range=attr(x, "year.range"))
  x$fact <- mkseas(x, width, start.day)
  days <- attr(x$fact, "bin.lengths")
  wet.ylim <- range(x$wet, na.rm=TRUE)
  dry.ylim <- range(x$dry, na.rm=TRUE)
  if (ylog)
    ylog <- "y"
  else {
    ylog <- ""
    wet.ylim[1] <- 0
    dry.ylim[1] <- 0
  }
  if (!missing(maxy)) {
    if (length(maxy) == 1)
      maxy[2] <- maxy
    wet.ylim[2] <- maxy[1]
    dry.ylim[2] <- maxy[2]
  }
  num.fact <- length(levels(x$fact))
  num <- num.fact + rep
  seas.boxplot <- function(x, at, var, col) {
    pl <- boxplot(as.formula(paste(var, "~ fact")), x, at=at, col=col, log=ylog,
                  varwidth=TRUE, names=NA, add=TRUE)
    lines(at, tapply(x[[var]], x$fact, mean, na.rm=TRUE), lwd=2)
    invisible(pl)
  }
  wet <- x
  dry <- x
  mar.wet <- c(2, 4, 4, 2) + 0.1
  mar.dry <- c(5, 4, 1, 2) + 0.1
  op <- getOption("seas.interarrival")
  par(mfrow=c(2, 1))
  frame()
  par(mar=mar.wet, yaxs="i", xaxs="r", bty="l")
  plot.window(c(0.5, num + 0.5), ylim=wet.ylim, log=ylog)
  .seasmonthgrid(width, days, start, rep, start.day)
  seas.boxplot(wet, 1:num.fact, "wet", op$col[1])
  if (rep > 0) {
    n <- 1
    if (rep >= num.fact) {  # whole repetitions
      for (n in 1:floor(rep / num.fact))
        seas.boxplot(wet, (num.fact * n + 1):(num.fact * (n + 1)),
                     "wet", op$col[1])
      n <- n + 1
    }
    if (rep %% num.fact != 0) {  # non-whole repetitions
      lv <- levels(wet$fact)[1:(rep %% num.fact)]
      wet <- wet[wet$fact %in% lv,]
      wet$fact <- factor(wet$fact, lv)
      seas.boxplot(wet, (num.fact * n + 1):(num.fact * n + rep %% num.fact),
                   "wet", op$col[1])
    }
  }
  axis(1, 1:num, ((1:num) + start - 2) %% num.fact + 1)
  title(main, xlab=xlab, ylab=ylab1)

  frame()
  par(mar=mar.dry, yaxs="i", xaxs="r", bty="c")
  plot.window(c(0.5, num + 0.5), ylim=rev(dry.ylim), log=ylog)
  .seasmonthgrid(width, days, start, rep, start.day, FALSE)
  seas.boxplot(dry, 1:num.fact, "dry", op$col[2])
  if (rep > 0) {
    n <- 1
    if (rep >= num.fact) {  # whole repetitions
      for (n in 1:floor(rep / num.fact))
        seas.boxplot(dry, (num.fact * n + 1):(num.fact * (n + 1)),
                     "dry", op$col[2])
      n <- n + 1
    }
    if (rep %% num.fact != 0) {  # non-whole repetitions
      lv <- levels(dry$fact)[1:(rep %% num.fact)]
      dry <- dry[dry$fact %in% lv,]
      dry$fact <- factor(dry$fact, lv)
      seas.boxplot(dry, (num.fact * n + 1):(num.fact * n + rep %% num.fact),
                   "dry", op$col[2])
    }
  }
  title(xlab=xlab, ylab=ylab2)
  axis(3, at=1:num, labels=FALSE)
}
