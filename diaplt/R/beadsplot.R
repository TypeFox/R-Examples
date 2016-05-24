# beads plot
# http://code.google.com/p/cowares-excel-hello/source/browse/trunk/beadsplot/
#
# Copyright (C) 2013 Tomizono
# Fortitudinous, Free, Fair, http://cowares.nobody.jp
#                            http://paidforeveryone.wordpress.com/

# beadsplot : S3 Method
#  default = data.frame
#  numeric
#  character
#  formula

# beadsplot.data.frame(
#   x = data.frame, observations
#   index = character vector, factor to separate rows, 
#           default=NULL is without factors
#   horizontal = logical, default=FALSE,
#                TRUE is to exchange x and y axes
#   col = colors to distinguish factors,
#         default=NULL is automatic coloring
#   sheer = alpha level to sheer colors,
#           default=NULL is c(1, 0.9, 0.6, 0.2)
#   shading = shading density to draw inside of beads,
#             default=NA is automatic, usually no shadings
#   shading.angle = shading angle to draw inside of beads,
#                   default=NA is automatic, usually no shadings
#   bw = numeric, beads width relative to series distance,
#        default=0.2
#   lwd = numeric vector, line width of beads border,
#         default=1
#   lwd.center =  numeric vector, line width of beads center,
#                 default=lwd, same width as border 
#   legend = logical, default=TRUE, to draw legend
#   label.factor = logical, default=TRUE, to draw labels of factor
#                  character vector, as alternative labels
#   label.range = logical, default=TRUE, to draw range value of series
#   drift.label.factor = numeric vector, (value, cycle),
#                        give small drifts on factor label location,
#                        default=c(0.2,2)
#   drift.label.range = numeric vector, (value, cycle),
#                       give small drifts on range label location,
#                       default=c(0,0), no drifts
#   S = summary function for bottom vertex,
#       default=min
#   E = summary function for center diagonal,
#       default=mean
#   N = summary function for top vertex,
#       default=max
#   summary.labels = character vector, c('min','mean','max')
#                    explicit names for summary functions,
#                    default=NULL, no explicit names
#   plot = logical, default=TRUE, to draw a chart
#   verbose = logical, default=FALSE, TRUE is to show debug information
#   ... = accepts graphical parameters and scale parameters

# beadsplot.numeric(
#   x = numeric vector, single observations

# beadsplot.character(
#   x = character, column name of factor
#   data = data.frame, contains a factor and observations

# beadsplot.formula(
#   formula = formula, a model formula, ex. factor~obs1+obs2+obs3
#   data = data.frame, contains variables in formula

# scale parameters
#
#   scale.range=1,          # width between center and border grids
#   scale.mean=0,           # location of center grid, NULL turns off
#   scale.log=F,            # logical value to enable log10 scaling
#   scale.data.center=NULL, # give center value vector from outside
#   scale.data.border=NULL, # give border value matrix from outside
#   scale.grid.center=NA,   # color of center grid, NULL turns off
#   scale.grid.border=NA,   # color of border grids, NULL turns off
#   cex.axis=1              # font size of grid label

# data structure (output values)
#
# list of 3 item (scaled, raw, scale)
#
# scaled and raw are same structure,
# raw is summary stats,
# scaled is conversion of raw, using scale
#
# three-dimensional array (series, factors, summaries)
#   series: each column of data
#   factors: each index
#   summaries: S/min, E/mean, N/max
#
# for one-dimensional vector data, series=1 is used.
# for one-level index=NULL, factors=1 is used.


beadsplot <- function(x, ...) UseMethod('beadsplot')

beadsplot.data.frame <- function(x, index=NULL, horizontal=FALSE,
               col=NULL, sheer=NULL, shading=NA, shading.angle=NA,
               bw=0.2, lwd=1, lwd.center=lwd, 
               legend=TRUE, label.factor=TRUE, label.range=TRUE, 
               drift.label.factor=c(0.2,2), drift.label.range=c(0,0),
               S=min, E=mean, N=max, summary.labels=NULL,
               plot=TRUE, verbose=FALSE, ...) {

  if(is.null(index)) {
    index <- rep('', nrow(x))
  } else if(!is.null(dim(index))) {
    index <- as.character(as.vector(simplify2array(index)))
  } else {
    index <- as.character(index)
  }

  scale <- make.scale(...)
  stats <- calc.beads(x, index, scale, S, E, N, summary.labels, verbose)
  #attr(stats, 'call') <- match.call(expand.dots=TRUE)
  if(verbose) print(stats)
 
  if(plot) {
    plt.beads(stats=stats, scale=scale, 
              horizontal=horizontal, bw=bw, col=col, 
              shading=shading, shading.angle=shading.angle, 
              lwd=lwd, lwd.center=lwd.center, legend=legend, 
              label.factor=label.factor, label.range=label.range, 
              drift.label.factor=drift.label.factor, 
              drift.label.range=drift.label.range, 
              sheer=sheer, summary.labels=summary.labels,
              verbose=verbose, ...)
    invisible(stats)
  } else {
    stats
  }
}

beadsplot.numeric <- function(x, index=NULL, horizontal=FALSE,
               col=NULL, sheer=NULL, shading=NA, shading.angle=NA,
               bw=0.2, lwd=1, lwd.center=lwd, 
               legend=TRUE, label.factor=TRUE, label.range=TRUE, 
               drift.label.factor=c(0.2,2), drift.label.range=c(0,0),
               S=min, E=mean, N=max, summary.labels=NULL,
               plot=TRUE, verbose=FALSE, ...) {
  beadsplot.data.frame(data.frame(x),
    index=index, horizontal=horizontal, col=col, sheer=sheer, 
    shading=shading, shading.angle=shading.angle,
    bw=bw, lwd=lwd, lwd.center=lwd.center, legend=legend, 
    label.factor=label.factor, label.range=label.range,
    drift.label.factor=drift.label.factor, 
    drift.label.range=drift.label.range, 
    S=S, E=E, N=N, summary.labels=summary.labels,
    plot=plot, verbose=verbose, ...)
}

beadsplot.character <- function(x, data, horizontal=FALSE,
               col=NULL, sheer=NULL, shading=NA, shading.angle=NA,
               bw=0.2, lwd=1, lwd.center=lwd, 
               legend=TRUE, label.factor=TRUE, label.range=TRUE, 
               drift.label.factor=c(0.2,2), drift.label.range=c(0,0),
               S=min, E=mean, N=max, summary.labels=NULL,
               plot=TRUE, verbose=FALSE, ...) {
  # convert column name to column number to use - operator
  column <- which(names(data) == x)
  if(length(column) == 0) column <- as.numeric(x)
  
  beadsplot.data.frame(data[- column], index=data[,column], 
    horizontal=horizontal, col=col, sheer=sheer, 
    shading=shading, shading.angle=shading.angle,
    bw=bw, lwd=lwd, lwd.center=lwd.center, legend=legend, 
    label.factor=label.factor, label.range=label.range,
    drift.label.factor=drift.label.factor, 
    drift.label.range=drift.label.range, 
    S=S, E=E, N=N, summary.labels=summary.labels,
    plot=plot, verbose=verbose, ...)
}

beadsplot.formula <- function(formula, data, horizontal=FALSE,
               col=NULL, sheer=NULL, shading=NA, shading.angle=NA,
               bw=0.2, lwd=1, lwd.center=lwd, 
               legend=TRUE, label.factor=TRUE, label.range=TRUE, 
               drift.label.factor=c(0.2,2), drift.label.range=c(0,0),
               S=min, E=mean, N=max, summary.labels=NULL,
               plot=TRUE, verbose=FALSE, ...) {
  # generate a model data frame,
  # with the 1st column as factor, and the 2nd- as observations
  x <- model.frame(formula, data=data)
  
  beadsplot.data.frame(x[-1], index=x[1], 
    horizontal=horizontal, col=col, sheer=sheer, 
    shading=shading, shading.angle=shading.angle,
    bw=bw, lwd=lwd, lwd.center=lwd.center, legend=legend, 
    label.factor=label.factor, label.range=label.range,
    drift.label.factor=drift.label.factor, 
    drift.label.range=drift.label.range, 
    S=S, E=E, N=N, summary.labels=summary.labels,
    plot=plot, verbose=verbose, ...)
}

beadsplot.default <- beadsplot.data.frame


blank.stats <- function(x, index) {
  # x as data.frame, index as vector of character
  series <- names(x)
  factors <- sort(unique(index))
  summaries <- c('S', 'E', 'N')
  # (bottom vertex, center diagonal, top vertex)
  array(NA, 
        dim=c(series=length(series), 
              factors=length(factors), 
              summaries=length(summaries)),
        dimnames=list(series=series, 
                      factors=factors, 
                      summaries=summaries))
}

calc.stats <- function(stats, x, index, S, E, N, 
                       summary.labels=NULL, verbose=FALSE) {
  # stats as a blank.stats
  # x as data.frame, index as vector of character

  SUMMARY <- list(S, E, N)
  names(SUMMARY) <- dimnames(stats)$summaries
  if(verbose) str(SUMMARY)

  if(is.null(summary.labels))
    summary.labels <- sapply(SUMMARY, 
                             function(a) paste(deparse(a), collapse=' '))
  if(verbose) print(summary.labels)
  attr(stats, 'summary.labels') <- summary.labels

  nseries <- dim(stats)['series']
  nsummaries <- dim(stats)['summaries']

  for(series in 1L:nseries) {
    for(summaries in 1L:nsummaries) {
      stats[series,,summaries] <- 
        tapply(x[,series], index, 
               function(a) SUMMARY[[summaries]](a[!is.na(a)])
              )
    }
  }

  stats
}

scale.stats.range <- function(stats, zoom=2, mid=0, centers=NULL) {
  lseries <- as.list(1L:dim(stats)['series'])
  whole.range <- sapply(lseries, 
                        function(series) range(stats[series,,c('S','N')]))
  if(is.null(mid)) {
    mid <- 0
    whole.mean <- apply(whole.range, 2, mean)
  } else if(is.null(centers)) {
    whole.mean <- sapply(lseries, 
                         function(series) mean(stats[series,,'E']))
  } else {
    whole.mean <- centers
  }
  scale.stats.any(stats, zoom, mid, whole.range, whole.mean)
}

scale.stats.any <- function(stats, zoom, mid, intervals, centers=NULL) {
  # intervals must be matrix or data.frame 
  #  with 2 rows (min, max) and series columns
  # centers must be a vector for series
  nseries <- dim(stats)['series']
  for(series in 1L:nseries) {
    center <- if(is.null(centers)) mean(intervals[,series]) else centers[series]
    scale <- diff(intervals[,series]) / zoom
    if(!is.finite(scale)) scale <- 1
    stats[series,,] <- (stats[series,,] - center) / scale + mid
  }
  stats
}

scale.stats <- function(stats, x, index, scale=make.scale()) {
  mid <- scale$scale.mean
  zoom <- scale$scale.range
  if(!is.null(zoom)) zoom <- 2 * zoom

  if(is.null(scale$scale.data.border)) {
    if(!is.null(zoom))
      stats <- scale.stats.range(stats, zoom, mid,
                                 scale$scale.data.center)
  } else {
    if(is.null(zoom)) zoom <- 2
    if(is.null(mid)) mid <- 0
    stats <- scale.stats.any(stats, zoom, mid, 
                             scale$scale.data.border, 
                             scale$scale.data.center)
  }

  if(scale$scale.log) {
    saturate.log <- 0.1
    stats[stats < saturate.log] <- saturate.log
    stats <- log10(stats)
  }

  stats
}

calc.beads <- function(x, index, scale, S=min, E=mean, N=max, 
                       summary.labels=NULL, verbose=FALSE) {
  # x as data.frame, index as vector of character
  stats <- blank.stats(x, index)
  stats <- calc.stats(stats, x, index, S, E, N, summary.labels, verbose)
  scaled <- scale.stats(stats, x, index, scale)
  list(scaled=scaled, raw=stats, scale=scale)
}

plt.gridline <- function(scale, xy) {
  labels <- c(C='C', U='U', L='L')

  pt <- matrix(nrow=length(labels), ncol=2, dimnames=list(labels, xy))
  side <- rep(NA, 2)
  names(side) <- xy
  side[c('x','y')] <- c(4, 1)

  center <- scale$scale.mean
  border <- scale$scale.range
  col.center <- scale$scale.grid.center 
  col.border <- scale$scale.grid.border 

  draw.grids <- c(center=!is.null(center), border=!is.null(border)) &
                c(!is.null(col.center), !is.null(col.border))

  grids.at <- c(C=if(draw.grids['center']) center else NA)
  grids.at[c('U','L')] <- if(draw.grids['border']) {
                            c(1, -1) * border +
                              if(is.null(center)) 0 else center
                          } else {
                            rep(NA, 2)
                          }
  if(scale$scale.log) grids.at <- log10(grids.at)

  grids.col <- c(C=col.center, U=col.border, L=col.border)

  pt[,'x'] <- NA
  pt[,'y'] <- grids.at

  for(cul in names(labels)) {
    if(is.finite(grids.at[cul])) { # reject +-Inf and NA
      abline(v=pt[cul,1], h=pt[cul,2], col=grids.col[cul])
      mtext(text=labels[cul], side=side[1], 
            at=grids.at[cul], col=grids.col[cul], 
            cex=scale$cex.axis, las=1)
    }
  }
}

make.scale <- function(...) {
  # extract scale parameters from dots args

  # paramter defaults
  scale.init <- 
    list(scale.range=1,          # width between center and border grids
         scale.mean=0,           # location of center grid, NULL turns off
         scale.log=F,            # logical value to enable log10 scaling
         scale.data.center=NULL, # give center value vector from outside
         scale.data.border=NULL, # give border value matrix from outside
         scale.grid.center=NA,   # color of center grid, NULL turns off
         scale.grid.border=NA,   # color of border grids, NULL turns off
         cex.axis=1              # font size of grid label
        )
  color.default <- 'grey'

  scale <- modifyList(scale.init, list(...))

  if(!is.null(scale$scale.grid.center) && is.na(scale$scale.grid.center)) {
    scale$scale.grid.center <- 
      if(is.null(scale$scale.mean)) NULL else color.default
  }
  if(!is.null(scale$scale.grid.border) && is.na(scale$scale.grid.border)) { 
    scale$scale.grid.border <- 
      if(is.null(scale$scale.range)) NULL else color.default
  }

  scale
}

plt.init <- function(lim, ...) {
  axt <- rep('s', 2)
  names(axt) <- dimnames(lim)[[2]]
  axt['x'] <- 'n'
  do.call.with.par('plot.default', list(...), 
                   xlim=lim[,1], ylim=lim[,2], 
                   x=NA, xlab='', ylab='', xaxt=axt[1], yaxt=axt[2])
}

calc.drifts <- function(scaled, series, 
                        label.factor, drift.label.factor,
                        label.range, drift.label.range) {
  # generate drift values for label position
  drifts <- list(f=NULL, r=NULL)

  if(label.factor) {
    rank.f <- order(order(scaled[series,,'E']))
    drifts$f <- calc.drifts.vector(rank.f, drift.label.factor)
  }
  if(label.range) {
    drifts$r <- calc.drifts.vector(series, drift.label.range)
  }

  drifts
}

calc.drifts.vector <- function(data, seed) {
  if(is.null(seed)) seed <- NA
  s1 <- seed[1]
  s2 <- seed[2]
  drift <-
    if(is.na(s2)) {
      rep(s1, length(data))
    } else {
      data %% s2 * s1
    }
  drift[is.nan(drift) | is.na(drift)] <- 0
  drift
}

make.colors <- function(n, col=NULL, sheer=NULL) {
  # generate color matrix with n columns and sheer rows
  if(is.null(sheer)) sheer=c(1, 0.9, 0.6, 0.2)
  nrow <- length(sheer)
  col <- if(is.null(col)) rainbow(n) else 
           if(length(as.vector(col)) != n * nrow)
             rep(col, n)[1L:n] else
               col
  
  cols <- matrix(data=col, ncol=n, nrow=nrow, byrow=T)
  sheers <- matrix(sheer, ncol=n, nrow=nrow)
  colors <- matrix(ncol=n, nrow=nrow)

  for(factors in 1L:n) {
    colors[,factors] <- 
      sapply(as.list(1L:nrow), function(r) 
        adjustcolor(cols[r,factors], alpha.f=sheers[r,factors]))
  }

  colors
}

plt.label.factors <- function(text, pt, drift,
                              col.label,
                              ...) {
  # draw factor labels at each diamond
  pt[,'x'] <- pt[,'x'] + drift
  do.call.with.par('text', list(...), dots.win=T,
                   x=pt[1,1], y=pt[3,2], col=col.label,
                   labels=text)
}

plt.label.range <- function(series, xy, drift, ra, sc, ...) {
  # draw labels at range of each series
  pt <- matrix(nrow=2, ncol=2, dimnames=list(c('min','max'), xy))
  pt[,'x'] <- series
  pt['min','y'] <- min(sc[series,,'S']) - drift 
  pt['max','y'] <- max(sc[series,,'N']) + drift 
  do.call.with.par('text', list(...), dots.win=T,
                   x=pt['min',1], y=pt['min',2], 
                   labels=min(ra[series,,'S']))
  do.call.with.par('text', list(...), dots.win=T,
                   x=pt['max',1], y=pt['max',2], 
                   labels=max(ra[series,,'N']))
}

plt.label.series <- function(text, horizontal, ...) {
  # draw labels on series axis
  side <- if(horizontal) 2 else 1
  at <- 1L:length(text)
  do.call.with.par('axis', list(...), side=side, at=at, labels=text)
}

plt.legend <- function(text, col, lim, ...) {
  # draw factor names above the chart
  nfactors <- length(text)
  pt <- seq(from=lim[1,1], to=lim[2,1], length.out=nfactors)
  for(factors in 1L:nfactors) {
    do.call.with.par('mtext', list(...), dots.win=T,
                     text=text[factors], 
                     side=3, at=pt[factors], col=col[factors])
  }
}

plt.legend.summary <- function(text, col='black', ...) {
  # explicit names for summary functions are shown
  if(!is.null(text))
    do.call.with.par('mtext', list(...), dots.win=T,
                     text=paste(text, collapse=', '), 
                     side=3, line=1, col=col)
}

plt.diamond <- function(pt, col, border, lwd=1, lwd.center=lwd,
                        density=NULL, angle=45) {
  # draw a rhombus with a diagonal of center
  #
  # pt is numerical matrix of 4 rows 2 columns.
  # 1st column is x, 2nd is y, always.  ignore names.
  #
  # x and y are vectors of 4 numerics, which means,
  #  center-right, top, center-left, bottom
  #
  # lwd can be NULL to turn off border lines
  # lwd.center can be NULL to turn off center diagonal

  polygon(x=pt[,1], y=pt[,2], col=col, lwd=lwd, 
          density=density, angle=angle,
          border=if(is.null(lwd)) NA else border)
  if(!is.null(lwd.center))
    lines(x=pt[c(1,3),1], y=pt[c(1,3),2], col=border, lwd=lwd.center)
}

plt.beads <- function(stats, scale, 
                      horizontal=FALSE,
                      bw=0.2, col=NULL, 
                      shading=NA, shading.angle=NA,
                      lwd=1, lwd.center=lwd,
                      legend=TRUE, label.factor=TRUE, label.range=TRUE, 
                      drift.label.factor=c(0.2,2),
                      drift.label.range=c(0,0),
                      sheer=NULL, summary.labels=NULL,
                      verbose=FALSE, ...) {
  xy <- c('x','y')
  if(horizontal) xy <- rev(xy)
  # x here means always series, y is always observations.
  # plt.* routines ignore the names, draw 1st as x-axis and 2nd as y-axis.
  # this trick support horizontal drawing switching x and y axes.

  pt <- matrix(nrow=4, ncol=2, dimnames=list(NULL, xy))
  lim <- matrix(nrow=2, ncol=2, dimnames=list(NULL, xy))

  sc <- stats$scaled
  ra <- stats$raw

  nseries <- dim(sc)['series']
  nfactors <- dim(sc)['factors']

  lw <- rep(lwd, nfactors)
  lwc <- rep(lwd.center, nfactors)

  colors <- make.colors(nfactors, col, sheer)
  dimnames(colors) <- list(c('legend','label','border','diamond'),
                           NULL)
  shadings <- make.shadings(nfactors, shading, shading.angle, verbose)

  # label.factor can be character vecotor for abbreviation
  if(is.character(label.factor) && length(label.factor) == nfactors) {
    factor.labels <- label.factor
    legend.labels <- paste(factor.labels, dimnames(sc)$factors, sep=': ')
  } else {
    factor.labels <- dimnames(sc)$factors
    legend.labels <- factor.labels
  }
  if(!is.logical(label.factor)) label.factor=TRUE

  lim[,'x'] <- c(1 - bw, nseries + bw)
  lim[,'y'] <- range(as.vector(sc))

  plt.init(lim=lim, ...)

  plt.gridline(scale=scale, xy=xy)

  for(series in 1L:nseries) {
    drifts <- calc.drifts(sc, series, label.factor, drift.label.factor,
                          label.range, drift.label.range)

    for(factors in 1L:nfactors) {
      # x and y of pt are vectors of 4 numerics, which means,
      #  center-right, top, center-left, bottom
      pt[,'x'] <- c(bw, 0, -bw, 0) + series
      pt[,'y'] <- sc[series,factors,][c('E','N','E','S')]

      plt.diamond(pt=pt, 
                  col=colors['diamond',factors], 
                  border=colors['border',factors], 
                  lwd=lw[factors], lwd.center=lwc[factors],
                  density=shadings$density[factors], 
                  angle=shadings$angle[factors])

      if(label.factor) 
        plt.label.factors(text=factor.labels[factors],
                          pt=pt, drift=drifts$f[factors],
                          col.label=colors['label',factors], ...)
    }
    
    if(label.range) 
      plt.label.range(series=series, xy=xy, drift=drifts$r,
                      ra=ra, sc=sc, ...)
  }

  plt.label.series(text=dimnames(sc)$series, 
                   horizontal=horizontal, ...)
  
  if(legend) {
    plt.legend(text=legend.labels,
               col=colors['legend',], lim=lim, ...)
    plt.legend.summary(text=summary.labels, ...)
  }
}

