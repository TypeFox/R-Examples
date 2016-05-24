make.par.list <- function(nams, x, n, m, def, recycle = sum(unnamed) > 0) {
##FIXME: should defaults for n, m, def be available?

# if nams are the names of our variables and x is a parameter
# specification such as list(a = c(1,2), c(3,4)) then 
# create a new list which uses the named variables from x
# and assigns the unnamed in order.  For the remaining variables
# assign them the default value if recycle = FALSE or recycle the
# unnamed variables if recycle = TRUE.  The default value for
# recycle is TRUE if there is at least one unnamed variable
# in x and is false if there are only named variables in x.
# n is the length of the series and m is the total number of series
# It only needs to know whether m is 1 or greater than m.
# def is the default value used when recycle = FALSE
# recycle = TRUE means recycle unspecified values
# recycle = FALSE means replace values for unspecified series with def
# Within series recycling is done even if recycle=FALSE.
  # Should we allow arbirary names in 1d case?
  # if (m > 1) stopifnot(all(names(x) %in% c("", nams)))
  if (!is.list(x)) x <- if (m == 1) list(x) else as.list(x)
  y <- vector(mode = "list", length = length(nams))
  names(y) <- nams
  in.x <- nams %in% names(x)
  unnamed <- if (is.null(names(x))) rep(TRUE, length(x)) else names(x) == ""
  if (!recycle) y[] <- def
  y[in.x] <- x[nams[in.x]]
  if (recycle) {
    stopifnot(sum(unnamed) > 0)
    y[!in.x] <- rep(x[unnamed], length.out = sum(!in.x)) ## CHECK, this was: x[unnamed]
  } else {
    y[which(!in.x)[seq_len(sum(unnamed))]] <- x[unnamed]
  }
  lapply(y, function(y) if (length(y)==1) y else rep(y, length.out = n))
}

plot.zoo <- function(x, y = NULL, screens, plot.type, panel = lines, 
  xlab = "Index", ylab = NULL, main = NULL, xlim = NULL, ylim = NULL,
  xy.labels = FALSE, xy.lines = NULL, yax.flip = FALSE,
  oma = c(6, 0, 5, 0), mar = c(0, 5.1, 0, if(yax.flip) 5.1 else 2.1), 
  col = 1, lty = 1, lwd = 1, pch = 1, type = "l", log = "",
  nc, widths = 1, heights = 1, ...)
{
  ## if y supplied: scatter plot y ~ x
  if(!is.null(y)) {
    if(NCOL(x) > 1 || NCOL(y) > 1) stop("scatter plots only for univariate zoo series")
    xyzoo <- merge.zoo(x, y, all = FALSE)
    xy <- coredata(xyzoo)
    xy <- xy.coords(xy[,1], xy[,2])

    xlab <- if(missing(xlab)) deparse(substitute(x)) else xlab
    ylab <- if(missing(ylab)) deparse(substitute(y)) else ylab
    xlim <- if(is.null(xlim)) range(xy$x[is.finite(xy$x)]) else xlim
    ylim <- if(is.null(ylim)) range(xy$y[is.finite(xy$y)]) else ylim
    if(is.null(main)) main <- ""
    do.lab <- if(is.logical(xy.labels)) xy.labels else TRUE
    if(is.null(xy.lines)) xy.lines <- do.lab
    ptype <- if(do.lab) "n" else if(missing(type)) "p" else type

    plot.default(xy, type = ptype,col = col, pch = pch, main = main,
      xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, log = log, ...)
    if(do.lab) text(xy, col = col,
      labels = if(!is.logical(xy.labels)) xy.labels else index2char(index(xyzoo)), ...)
    if(xy.lines) lines(xy, col = col, lty = lty, lwd = lwd, type = if(do.lab) "c" else "l", ...)

    return(invisible(xyzoo))
  }
  ## Else : no y, only x

  recycle <- function(a, len, nser)
     rep(lapply(as.list(a), rep, length.out = len), length.out = nser)
  # same as range except it passes pairs through
  range2 <- function(x, ...) if (length(x) == 2) x else range(x, ...)
  if (missing(plot.type)) {
	plot.type <- if (missing(screens)) "multiple"
	else if (length(unique(screens) == 1)) "single" 
	else "multiple"
  }
  plot.type <- match.arg(plot.type, c("multiple", "single"))
  nser <- NCOL(x)
  if (missing(screens)) {
	screens <- if (plot.type == "single") 1 else seq_len(nser)
  }	
  dots <- list(...)
  x.index <- index(x)
  if(is.ts(x.index)) x.index <- as.vector(x.index)
  cn <- if (is.null(colnames(x))) paste("V", seq_len(nser), sep = "")
	  else colnames(x)

  screens <- make.par.list(cn, screens, NROW(x), nser, 1)
  screens <- as.factor(unlist(screens))[drop = TRUE]
  ngraph <- length(levels(screens))
  if(nser > 1 && (plot.type == "multiple" || ngraph > 1)) {
    if (ngraph == 1) { 
	screens <- as.factor(seq(nser))
	ngraph <- nser
    }
    if(is.null(main)) main <- deparse(substitute(x))
    main.outer <- TRUE
    if(is.null(ylab)) ylab <- colnames(x)[!duplicated(screens)]
    if(is.null(ylab)) ylab <- paste("Series", which(!duplicated(screens)))
    if(is.call(ylab)) ylab <- as.expression(ylab)
    ylab <- rep(ylab, length.out = ngraph)
    if(!is.list(ylab)) ylab <- as.list(ylab)
    lty <- rep(lty, length.out = nser)
    lwd <- rep(lwd, length.out = nser)
    col <- make.par.list(cn, col, NROW(x), nser, 1)
    pch <- make.par.list(cn, pch, NROW(x), nser, par("pch"))
    type <- make.par.list(cn, type, NROW(x), nser, "l")
    if (!is.null(ylim)) {
        if (is.list(ylim)) ylim <- lapply(ylim, range2, na.rm = TRUE)
	else ylim <- list(range2(ylim, na.rm = TRUE))
	ylim <- lapply(make.par.list(cn, ylim, 2, nser, NULL), function(x) 
		if (is.null(x) || length(na.omit(x)) == 0) NULL 
		else range2(x, na.rm = TRUE))
    }
    panel <- match.fun(panel)
    if(missing(nc)) nc <- if(ngraph >  4) 2 else 1
    oldpar <- par(no.readonly = TRUE)
    on.exit({ par(oldpar) })
    nr <- ceiling(ngraph / nc)
    layout(matrix(seq(nr*nc), nr), widths = widths, heights = heights)
    par(mar = mar, oma = oma)
	# TRUE if all elements of L are the same -- else FALSE
	allsame <- function(L) {
		f <- function(x, y) if (identical(x, y)) x
		!is.null(Reduce(f, L))
	}
	# idx is vector of indices into ylim.  
	# If the entries indexed by it are all the same then use that common value;
	# otherwise, if the ylim are specified use the range of the ylim values;
	# otherwise, use the range of the data
	f <- function(idx) if (allsame(ylim)) ylim[idx][[1]]
		else if (!is.null(ylim) && length(idx) > 0 && 
			length(unlist(ylim[idx])) > 0) range(ylim[idx], finite = TRUE)
		else range(x[, idx], na.rm = TRUE)
	# ranges is indexed by screen
	ranges <- tapply(1:ncol(x), screens, f)
    for(j in seq_along(levels(screens))) {
      panel.number <- j
      y.side <- if (j %% 2 || !yax.flip) 2 else 4
      range. <- rep(ranges[[j]], length.out = length(time(x)))
      if(j%%nr==0 || j == length(levels(screens))) {
			args <- list(x.index, range., xlab = "", ylab = "", yaxt = "n",
				xlim = xlim, ylim = ylim[[j]], log = log, ...)
			args$type <- "n"
			do.call("plot", args)
			mtext(xlab, side = 1, line = 3)
      } else {      
			args <- list(x.index, range., xaxt = "n", yaxt = "n", xlab = "", 
				ylab = "", xlim = xlim, ylim = ylim[[j]], log = log, ...)
			args$type <- "n"
			do.call("plot", args)
			if ("bty" %in% names(args) && args$bty == "n") {} else box()
      }
      do.call("axis", c(list(side = y.side, xpd = NA), dots))
      mtext(ylab[[j]], y.side, line = 3)

      for(i in which(screens == levels(screens)[j])) {
        ## for potential usage in panel function
        series.number <- i
        series.within.screen <- ave(seq_along(screens), screens, FUN = seq_along)[series.number]
       
        ## draw individual lines/points with panel function
        panel(x.index, x[, i], col = col[[i]], pch = pch[[i]], lty = lty[i], lwd = lwd[i], type = type[[i]], ...)
      }
    }
  } else {
    if(is.null(ylab)) ylab <- deparse(substitute(x))
    if(is.call(ylab)) ylab <- as.expression(ylab)
    if(is.null(main)) main <- ""
    main.outer <- FALSE
    if(is.null(ylim)) ylim <- range(x, na.rm = TRUE)
	else ylim <- range2(c(ylim, recursive = TRUE), na.rm = TRUE)

    lty <- rep(lty, length.out = nser)
	lwd <- rep(lwd, length.out = nser)
    col <- make.par.list(cn, col, NROW(x), nser, 1)
    pch <- make.par.list(cn, pch, NROW(x), nser, par("pch"))
    type <- make.par.list(cn, type, NROW(x), nser, "l")
   
    dummy <- rep(range(x, na.rm = TRUE), 
	length.out = length(index(x)))
	    
    args <- list(x.index, dummy, xlab = xlab, ylab = ylab[1], ylim = ylim, xlim = xlim, log = log, ...)
    args$type <- "n"
    do.call("plot", args)
	if ("bty" %in% names(args) && args$bty == "n") {} else box()
    y <- as.matrix(x)
    for(i in 1:nser) {
      panel(x.index, y[, i], col = col[[i]], pch = pch[[i]], lty = lty[i], 
        lwd = lwd[i], type = type[[i]], ...)
    }
  }
  dots <- list(...)
  title.args <- c(list(main = main, outer = main.outer),
    dots[grep("[.]main$", names(dots))])
  do.call("title", title.args)
  return(invisible(x))
}

lines.zoo <- function(x, y = NULL, type = "l", ...)
{
  if (is.null(y)) {
     if(NCOL(x) == 1) lines(index(x), x, type = type, ...)
       else stop("Can't plot lines for multivariate zoo object")
  } else
     lines(coredata(cbind(x,y)), type = type, ...)
}

points.zoo <- function(x, y = NULL, type = "p", ...)
  lines(x, y, type = type, ...)


plot.tis <- function(x, ...) eval.parent(substitute(plot(as.zoo(x), ...)))

plot.ti <- function (x, y, xlab = "", ...) 
{
	x <- tis::POSIXct(x)
	NextMethod()
}

points.ti <- lines.ti <- function(x, ...) {
	x <- tis::POSIXct(x)
	NextMethod()
}

