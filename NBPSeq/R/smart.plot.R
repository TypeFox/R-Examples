##' An alternative to plot.default() for plotting a large number of
##' densely distributed points.  This function can produce a visually
##' almost identical plot using only a subset of the points.  This is
##' particular useful for reducing output file size when plots are
##' written to eps files.
##'
##'   Writing plots with a large number of points to eps files can
##'   result in big files and lead to very slow rendering time.
##'
##'   Usually for a large number of points, a lot of them will
##'   overlap with each other. Plotting only a subset of selected
##'   non-overlapping points can give visually almost identical
##'   plots. Further more, the plots can be enhanced if using gray
##'   levels (the default setting) that are proportional to the
##'   number points overlapping with each plotted point.
##' 
##'   This function scans the points sequentially. For each unmarked
##'   point that will be plotted, all points that overlap with it
##'   will be marked and not to plotted, and the number of
##'   overlapping points will be recorded. This is essentially
##'   producing a 2d histogram. The freqs of the points will be
##'   converted to gray levels, darker colors correspond to higher
##'   freqs.
##'
##' @title (private) An alternative to plot.default() for plotting a large number of
##' densely distributed points.
##'
##' @param x  x
##' @param y  y
##' @param xlim xlim
##' @param ylim ylim
##' @param xlab x label
##' @param ylab y label
##' @param log log
##' @param resolution a number, determines the distance below which
##' points will be considered as overlapping.
##' @param plot logical, whether
##' @param col color
##' @param clip clip
##' @param color.clipped color of clipped points
##' @param ... other arguments are the same as in plot.default().
##'
##'
##' @return (if plot=FALSE) a list 
##'
##'   \item{x, y}{the x, y-coordinates of the subset of representative points}
##' 
##'   \item{id}{the indicies of these points in the original data set}
##'  
##'   \item{freqs}{the numbers of points that overlap with each representative point}
##'
##'   \item{col}{colors determined by the freqs}
smart.plot = function(x, y=NULL, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
  log="", resolution=100, plot=TRUE, col=NULL, clip=Inf, color.clipped=TRUE, ...) {


  ## These lines are copied from plot.default.
  xlabel = if (!missing(x)) deparse(substitute(x));
  ylabel = if (!missing(y)) deparse(substitute(y));
  xy = xy.coords(x, y, xlabel, ylabel, log);
  xlab = if (is.null(xlab)) xy$xlab else xlab;
  ylab = if (is.null(ylab)) xy$ylab else ylab;
  xlim = if (is.null(xlim)) range(xy$x[is.finite(xy$x)]) else xlim;
  ylim = if (is.null(ylim)) range(xy$y[is.finite(xy$y)]) else ylim;
  
  x = xy$x;
  y = xy$y;
  n = length(x);
  id = is.finite(xy$x) & is.finite(xy$y);
  id[id & (x < xlim[1] | x > xlim[2] | y < ylim[1] | y > ylim[2])]=FALSE;

  logxy = strsplit(log, NULL)[[1L]];
  if ("x" %in% logxy) {
    x[id] = log(x[id]);
    epsx = diff(log(xlim)) / resolution;
  } else
    epsx = diff(xlim) / resolution;

  if ("y" %in% logxy) {
    y[id] = log(y[id]);
    epsy = diff(log(ylim)) / resolution; 
  } else
    epsy = diff(ylim) / resolution; 

  counts = rep(0, n);
  i = 1;

  ## Scan the points and select non-overlapping ones
  while (i < n) {
    if (id[i]) {
      ids = ((i+1):n)[id[(i+1):n]];
      overlap = ids[abs(x[ids] - x[i]) < epsx & abs(y[ids] - y[i]) < epsy];
      id[overlap] = FALSE;
      counts[i] = length(overlap) + 1;
    }
    i = i + 1;
  }

  ## Sort the data so that points representing more points will be plotted later.
  id = (1:n)[id];
  id = id[order(counts[id])];
  counts = counts[id];

  if (is.null(col)) {
    ## Convert counts of overlapping points to gray levels
    counts.clipped = counts;
    id.clipped = counts > clip;
    counts.clipped[id.clipped] = clip;
    col = gray((1 - counts.clipped / max(counts.clipped)) * 0.8);
    if (color.clipped) {
      col[id.clipped] = rgb(counts[id.clipped]/max(counts), 0, 0);
    }
  }

  if (plot) {
    plot(x=xy$x[id], y=xy$y[id], xlim=xlim, ylim=ylim, log=log, xlab=xlab, ylab=ylab, col=col, ...);
    invisible();
  }

  invisible(list(x=xy$x[id], y=xy$y[id], freqs=counts, col = col, id=id));
}


##' See description of \code{\link{smart.plot}} for more details.
##'
##' @title (private) An alternative to point.default() for plotting a large number of
##' densely distributed points.
##'
##' @param x x
##' @param y y
##' @param resolution a number, determines the distance below which
##' points will be considered as overlapping.
##' @param col color
##' @param clip clip
##' @param color.clipped color of clipped points
##' @param ... other arguments are the same as in plot.default().
##'
##' @return NULL
smart.points = function(x, y=NULL, 
  resolution=50, col=NULL, clip=Inf, color.clipped=TRUE, ...) {

  xy = xy.coords(x, y);

  xlim = range(xy$x[is.finite(xy$x)]);
  ylim = range(xy$y[is.finite(xy$y)]);

  x = xy$x;
  y = xy$y;
  n = length(x);
  id = is.finite(xy$x) & is.finite(xy$y);
  id[id & (x < xlim[1] | x > xlim[2] | y < ylim[1] | y > ylim[2])]=FALSE;

  log = "";
  
  logxy = strsplit(log, NULL)[[1L]];
  if ("x" %in% logxy) {
    x[id] = log(x[id]);
    epsx = diff(log(xlim)) / resolution;
  } else
    epsx = diff(xlim) / resolution;

  if ("y" %in% logxy) {
    y[id] = log(y[id]);
    epsy = diff(log(ylim)) / resolution; 
  } else
    epsy = diff(ylim) / resolution; 

  counts = rep(0, n);
  i = 1;

  ## Scan the points and select non-overlapping ones
  while (i < n) {
    if (id[i]) {
      ids = ((i+1):n)[id[(i+1):n]];
      overlap = ids[abs(x[ids] - x[i]) < epsx & abs(y[ids] - y[i]) < epsy];
      id[overlap] = FALSE;
      counts[i] = length(overlap) + 1;
    }
    i = i + 1;
  }

  ## Sort the data so that points representing more points will be plotted later.
  id = (1:n)[id];
  id = id[order(counts[id])];
  counts = counts[id];

  if (is.null(col)) {
    ## Convert counts of overlapping points to gray levels
    counts.clipped = counts;
    id.clipped = counts > clip;
    counts.clipped[id.clipped] = clip;
    col = gray((1 - counts.clipped / max(counts.clipped)) * 0.8);
    if (color.clipped) {
      col[id.clipped] = rgb(counts[id.clipped]/max(counts), 0, 0);
    }
  }

  points(x=xy$x[id], y=xy$y[id], col=col, ...);
  invisible();
}

example.smart.plot = function() {
  x = rnorm(20000);
  y = rnorm(20000);

  ## debug(smart.plot);
  ## debug(plot.default);
  ## undebug(plot.default);
  par(mfrow=c(2,3));
  plot(x, y, pch=19);
  ## Plot with translucent color
  plot(x, y, col=rgb(0, 0, 1, 0.1), pch=19);
  smart.plot(x, y, resolution=100, pch=19);

  smart.plot(x, y, pch=19);
  smart.plot(x, y, resolution = 100, pch="+");
  smart.plot(x, y, resolution = 100, pch="+", col=1);

  smart.plot(x, y, resolution=100, clip=5, pch=19);
  smart.plot(x, y, resolution= 100, pch=19, log="xy");

  obj = smart.plot(x, y, resolution= 100, pch=19);

  obj = smart.plot(x, y, resolution= 100, pch=19, plot=FALSE);


}
