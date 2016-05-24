###########################################################################/**
# @RdocFunction cnr
#
# @title "Simple creation of a CopyNumberRegions object"
#
# \description{
#  @get "title" containing a single region.
# }
#
# @synopsis
#
# \arguments{
#   \item{start, stop}{Two @numerics specifying the start and stop
#     location of the region.  If \code{stop} is @NULL, then
#     \code{start[2]} is used instead.}
#   \item{mean}{A @numeric specifying the mean level of the region.}
#   \item{chromosome}{An @integer specifying the chromosome ID.}
#   \item{xScale}{The default scaling parameter used for \code{start}
#     and \code{stop}, which is then also used for plotting genomic
#     locations.}
#   \item{...}{Additional arguments passed to
#     @see "aroma.core::CopyNumberRegions".}
# }
#
# \value{
#   Return an @see "aroma.core::CopyNumberRegions" object.
# }
#
# @examples "../incl/cnr.Rex"
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
cnr <- function(start, stop=NULL, mean=1, chromosome=1, xScale=1e6, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'start':
  start <- Arguments$getNumerics(start, length=c(1,2));

  # Argument 'start':
  if (is.null(stop)) {
    if (length(start) != 2) {
      throw("If argument 'stop' is not specified/NULL, then argument 'start' must be of length two: ", length(start));
    }
    stop <- start[2];
    start <- start[1];
  } else if (length(start) != 1) {
    throw("If argument 'stop' is specified, then argument 'start' must be a single position: ", length(start));
  }

  CopyNumberRegions(chromosome=chromosome,
                    start=xScale*start, stop=xScale*stop,
                    mean=mean, ...);
} # cnr()


setMethodS3("getXScale", "CopyNumberRegions", function(this, ...) {
  scale <- this$.xScale;
  if (is.null(scale))
    scale <- 1e-06;
  scale;
})

setMethodS3("getYScale", "CopyNumberRegions", function(this, ...) {
  scale <- this$.yScale;
  if (is.null(scale))
    scale <- 1;
  scale;
})



setMethodS3("plot", "CopyNumberRegions", function(x, joined=TRUE, xlab="Position", ylab="CN", ylim=c(0,5), xScale=getXScale(this), yScale=getYScale(this), ...) {
  # To please R CMD check
  this <- x;

  # Argument 'joined':
  joined <- Arguments$getLogical(joined);


  # Create an empty plot
  x <- xRange(this);
  y <- c(0,0);
  plot(xScale*x, yScale*y, type="n", ylim=ylim, xlab=xlab, ylab=ylab, ...);

  # Draw the levels
  if (joined) {
    lines(this, ...);
  } else {
    drawLevels(this, ...);
  }

  invisible(this);
})


##############################################################################
# HISTORY
# 2010-09-09
# o Added getXScale(), getYScale() and plot() for CopyNumberRegions.
# o Added cnr() to aroma.cn.
# 2010-07-19
# o Created.
##############################################################################
