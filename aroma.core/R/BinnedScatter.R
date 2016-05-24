###########################################################################/**
# @RdocClass BinnedScatter
#
# @title "The BinnedScatter class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{data}{A Nx2 @numaric @matrix.}
#   \item{density}{...}
#   \item{map}{...}
#   \item{params}{A @list of parameters.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @examples "../incl/BinnedScatter.Rex"
#
# @author
#
# \seealso{
#   The spatial density is estimated by internal functions of the
#   \pkg{smoothScatter} package.
# }
#*/###########################################################################
setConstructorS3("BinnedScatter", function(data=NULL, density=NULL, map=NULL, params=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(data)) {
    # Argument 'data':
    dim <- dim(data);
    if (dim[2] != 2) {
      throw("Argument 'data' must be a two-column matrix: ", dim[2]);
    }
  }


  extend(list(
    data=data,
    density=density,
    map=map,
    params=params
  ), "BinnedScatter");
})


setMethodS3("reorder", "BinnedScatter", function(x, orderBy="density", decreasing=FALSE, ...) {
  # To please R CMD check
  object <- x;

  o <- order(object[[orderBy]], decreasing=decreasing);
  object$data <- object$data[o,,drop=FALSE];
  object$density <- object$density[o];
  params <- object$params;
  params$orderBy <- orderBy;
  params$decreasing <- decreasing;
  object$params <- params;

  object;
})


setMethodS3("points", "BinnedScatter", function(x, ...) {
  # To please R CMD check
  object <- x;

  points(object$data, ...);
})


setMethodS3("plot", "BinnedScatter", function(x, ...) {
  # To please R CMD check
  object <- x;

  plot(object$data, ...);
})


setMethodS3("subset", "BinnedScatter", function(x, subset, ...) {
  # To please R CMD check
  object <- x;

  object$data <- object$data[subset,,drop=FALSE];
  object$density <- object$density[subset];
  object;
})


setMethodS3("subsample", "BinnedScatter", function(object, size=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  weightFcn <- function(object, ...) {
    w <- 1/object$density;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  n <- dim(object$data)[1];

  # Argument 'size':
  if (is.null(size)) {
    size <- n;
  } else {
    size <- Arguments$getNumeric(size, range=c(0, n));
    if (size < 1) {
      size <- round(size*n);
      if (size > n)
        size <- n;
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate sample weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  w <- weightFcn(object);

  # Standarize weights
  w <- w / sum(w, na.rm=TRUE);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Randomized sampling according to weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subset <- resample(1:length(w), size=size, prob=w, replace=FALSE);
  # Not needed anymore
  w <- NULL;
  res <- subset(object, subset=subset, ...);

  res;
}) # subsample()


setMethodS3("binScatter", "matrix", function(x, nbin=128, orderBy="density", decreasing=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'x':
  dim <- dim(x);
  if (dim[2] != 2) {
    throw("Argument 'x' must be a two-column matrix: ", dim[2]);
  }

  # Argument 'orderBy':
  if (!is.null(orderBy)) {
    orderBy <- match.arg(orderBy, c("density"));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify density estimator
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Used to be a workaround since it used to be in 'geneplotter'.
  ns <- getNamespace("grDevices");
  calcDensity <- get(".smoothScatterCalcDensity", mode="function", envir=ns);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate the (x,y) density
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate density only from finite data points
  ok <- which(is.finite(x[,1]) & is.finite(x[,2]));
  x <- x[ok,,drop=FALSE];
  # Not needed anymore
  ok <- NULL;
  map <- calcDensity(x, nbin=nbin);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Map each data point to a bin
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  xm <- map$x1;
  ym <- map$x2;

  nx <- length(xm);
  ny <- length(ym);
  dx <- x[,1]-xm[1];
  dy <- x[,2]-ym[1];
  w <- xm[nx] - xm[1];
  h <- ym[ny] - ym[1];
  ixm <- round(dx/w * (nx - 1));
  iym <- round(dy/h * (ny - 1));
  binIdx <- (1 + iym*nx + ixm);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the density at each data point
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dens <- map$fhat;
  idens <- dens[binIdx];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup return structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- BinnedScatter(
    data=x,
    density=idens,
    map=map,
    params=list(nbin=nbin)
  );

  # Reorder data?
  if (!is.null(orderBy)) {
    res <- reorder(res, orderBy=orderBy, decreasing=decreasing);
  }

  res;
}) # binScatter()



############################################################################
# HISTORY:
# 2010-11-04
# o ROBUSTNESS: Now subsample() for BinnedScatter utilizes resample().
# 2010-10-27
# o CLEANUP: Dropped outdated backup import to geneplotter in binScatter().
# 2009-05-09
# o UPDATED: Now binScatter() of BinnedScatter don't need the 'geneplotter'
#   package if R v2.9.0+ is used.
# 2008-11-26
# o Added Rdoc comments with an example.
# o Added constructor and reorder().
# 2008-11-14
# o Created.
############################################################################
