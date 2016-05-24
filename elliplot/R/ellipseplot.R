# ellipse plot
# http://code.google.com/p/cowares-excel-hello/source/browse/trunk/ellipseplot/
#
# Copyright (C) 2013 Tomizono
# Fortitudinous, Free, Fair, http://cowares.nobody.jp
#
# ellipseplot(
#             x= data frame for x-axis; factors and observations,
#             y= data frame for y-axis; factors and observations,
#             SUMMARY= function generating summaries to write contours,
#             SHEER= function adjusting color levels,
#             plot= TRUE for chart / FALSE for print summary,
#             verbose= TRUE for debugging information,
#             ...= accepts plot parameters
#            )
#
# other data types supported
#   numeric vector : no factors
#   matrix : same format as data.frame
#   list   : list of observation vectors, separated by factors
#
# support y=NULL
#   plot single axis data, separated by horizontal indexes
#
# requires midpoints

ellipseplot <- function(x, ...) UseMethod("ellipseplot")

ellipseplot.data.frame <- function(x, y=NULL, 
                   SUMMARY=ninenum, SHEER=sheer.color,
                   plot=TRUE, verbose=FALSE, ...) {
  if(plot) xaxt <- par('xaxt')
  if(is.null(y)) { # generate indexes for single axis plot
    xaxt <- 'n'
    y <- x
    f <- levels(as.factor(y[,1]))
    fl <- length(f)
    x <- data.frame(f=rep(f, 2), 
                    x=rep(1L:fl, 2) + 0.5 * rep(c(-1,1), each=fl))
    names(x)[2] <- names(y)[1]
  }

  stats <- calc.stats(x, y, SUMMARY)
  axes <- list(x=names(x), y=names(y))
  
  if(plot) {
    many.ellipses(stats, axes, SHEER, xaxt=xaxt, ...)
    invisible(stats)
  } else {
    stats
  }
}

ellipseplot.numeric <- function(x, y=NULL, 
                   SUMMARY=ninenum, SHEER=sheer.color,
                   plot=TRUE, verbose=FALSE, ...) {
  xd <- data.frame(o=rep('o',length(x)), x)
  yd <- if(is.null(y)) NULL else 
          data.frame(o=rep('o',length(y)), y)
  ellipseplot.data.frame(xd, yd, SUMMARY, SHEER, plot, verbose, ...)
}

ellipseplot.matrix <- function(x, y=NULL, 
                   SUMMARY=ninenum, SHEER=sheer.color,
                   plot=TRUE, verbose=FALSE, ...) {
  xd <- data.frame(x, stringsAsFactors=F)
  xd[,2] <- as.numeric(xd[,2])
  yd <- NULL
  if(!is.null(y)) {
    yd <- data.frame(y, stringsAsFactors=F)
    yd[,2] <- as.numeric(yd[,2])
  }
  ellipseplot.data.frame(xd, yd, SUMMARY, SHEER, plot, verbose, ...)
}

ellipseplot.list <- function(x, y=NULL, 
                   SUMMARY=ninenum, SHEER=sheer.color,
                   plot=TRUE, verbose=FALSE, ...) {
  sx <- 1L:length(x)
  nx <- names(x)
  if(is.null(nx)) nx <- as.character(sx)
  xf <- xv <- NULL
  for(i in sx) {
    xv <- c(xv, x[[i]])
    xf <- c(xf, rep(nx[i], length(x[[i]])))
  }
  xd <- data.frame(factor=xf, x=xv)

  yd <- NULL
  if(!is.null(y)) {
    sy <- 1L:length(y)
    ny <- names(y)
    if(is.null(ny)) ny <- as.character(sy)
    yf <- yv <- NULL
    for(i in sy) {
      yv <- c(yv, y[[i]])
      yf <- c(yf, rep(ny[i], length(y[[i]])))
    }
    yd <- data.frame(factor=yf, y=yv)
  }

  ellipseplot.data.frame(xd, yd, SUMMARY, SHEER, plot, verbose, ...)
}

ellipseplot.default <- ellipseplot.data.frame

# draw multiple ellipses of stat
# stats is a list of stat
many.ellipses <- function(stats, axes, SHEER=sheer.color, ...) {
  xy <- list('x', 'y')
  lims <- lapply(xy, function(a) 
                 c(min=min(sapply(stats, function(stat) min(stat[a]))), 
                   max=max(sapply(stats, function(stat) max(stat[a]))))
                 )
  names(lims) <- paste(xy, 'lim', sep='')

  lims$xlab <- axes$x[2]
  lims$ylab <- axes$y[2]

  pars <- modifyList(lims, c(list(...), x=NA))
  do.call('plot', pars)

  statnum <- length(stats)
  col <- ( if(hasArg(col)) list(...)$col
           else rainbow(statnum) )
  if(length(col) < statnum) col <- rep(col, statnum)
  name <- names(stats)
  if(is.null(name)) name <- as.character(1L:statnum)

  for(i in 1L:statnum) {
    ellipses(stats[[i]], name[i], col[i], SHEER)
  }
}

# draw ellipses of stat
# stat must have odd number of rows and ascending order
ellipses <- function(stat, name, col, SHEER=sheer.color) {
  xy <- list('x', 'y')
 
  boxes <- calc.abox(stat)
  boxnumber <- length(boxes)
  
  center <- lapply(xy,
                   function(a) boxes[[1]][a, 'center'])
  names(center) <- xy

  for(i in 1L:boxnumber) {
    pcol <- SHEER(col, i / boxnumber)
    anellipse(boxes[[i]], col=pcol, border=pcol)
  }

  text(center$x, center$y, name)
  mtext(name, side=3, at=center$x, col=col)
  mtext(name, side=4, at=center$y, col=col)
}

sheer.color <- function(col, level) {
  sheer <- level^2 * 0.5
  adjustcolor(col, alpha.f=sheer)
}

# draw a single inscribed ellipse to the specified box 
# accepts parameters for polygon()
anellipse <- function(abox, verbose=FALSE, ...) {
  axes <- rbind(rt=abox[,'high'] - abox[,'center'], 
                lb=abox[,'center'] - abox[,'low'])
  colnames(axes) <- c('x', 'y')
  if(verbose) print(axes)

  qx <- c('rt','lb','lb','rt')
  qy <- c('rt','rt','lb','lb')

  seed.ellipse <- data.frame(
    quadrant=1:4,
    startangle=seq(0, 2*pi, length=5)[1:4],
    endangle=seq(0, 2*pi, length=5)[2:5],
    xcenter=rep(abox['x', 'center'], 4),
    ycenter=rep(abox['y', 'center'], 4),
    xaxis=axes[qx, 'x'],
    yaxis=axes[qy, 'y']
  )
  if(verbose) print(seed.ellipse)

  x.ellipse <- as.vector(apply(seed.ellipse, 1, calc.ellipse.x))
  y.ellipse <- as.vector(apply(seed.ellipse, 1, calc.ellipse.y))
  if(verbose) { 
    str(x.ellipse)
    str(y.ellipse)
  }

  polygon(x.ellipse, y.ellipse, ...)
}

calc.ellipse <- function(center, axis, 
                         start=0, end=2*pi, length=100, 
                         FUNC=cos) {
  theta <- seq(start, end, length=length)
  axis * FUNC(theta) + center
}

calc.ellipse.x <- function(seed) {
  calc.ellipse(seed['xcenter'], seed['xaxis'], 
               seed['startangle'], seed['endangle'],
               FUNC=cos)
}

calc.ellipse.y <- function(seed) {
  calc.ellipse(seed['ycenter'], seed['yaxis'], 
               seed['startangle'], seed['endangle'],
               FUNC=sin)
}


# expect a data frame with 1st column factor and 2nd column data,
# for each x and y
calc.stats <- function(x, y, SUMMARY=ninenum, na.rm=TRUE) {
  factors <- sort(union(levels(as.factor(x[,1])), 
                        levels(as.factor(y[,1]))))
  stats <- lapply(as.list(factors), function(f) {
                  calc.stat(x[x[,1]==f,2],
                            y[y[,1]==f,2],
                            SUMMARY)
           })
  names(stats) <- factors
  calc.stats.na.rm(stats, na.rm)
}

calc.stats.na.rm <- function(x, na.rm) {
  if(na.rm) 
    for(f in names(x))
      if(any(apply(x[[f]], 2, function(a) all(is.na(a))))) 
        x[f] <- NULL
  x
}

calc.stat <- function(x, y, SUMMARY=ninenum) {
  data.frame(x=SUMMARY(x), y=SUMMARY(y))
}

calc.abox <- function(stat) {
  np <- nrow(stat) + 1
  center <- np / 2
  lapply(as.list(1L:(center - 1)), 
         function(i) {
           box <- rbind(stat[i,], stat[center,], stat[np - i,])
           rownames(box) <- c('low', 'center', 'high')
           t(box)
         })
}

