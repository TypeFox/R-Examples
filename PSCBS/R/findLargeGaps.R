###########################################################################/**
# @RdocDefault findLargeGaps
# @alias findLargeGaps.data.frame
#
# @title "Identifies gaps of a genome where there exist no observations"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{chromosome}{(Optional) An @integer @vector of length J of
#     chromosome indices.}
#   \item{x}{A @numeric @vector of J of genomic locations.}
#   \item{minLength}{A positive @numeric scalar specifying the minimum
#     length of a gap.}
#   \item{resolution}{A non-negative @numeric specifying the minimum
#     length unit, which by default equals one nucleotide/base pair.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns @data.frame zero or more rows and with columns
#   \code{chromosome} (if given), \code{start}, \code{stop},
#   and \code{length}.
# }
#
# @author "HB"
#
# \seealso{
#   Use @see "gapsToSegments" to turn the set of identified gaps into
#   the complementary set of segments such that they can be passed
#   to @see "segmentByCBS", @see "segmentByPairedPSCBS" and
#   @see "segmentByNonPairedPSCBS" via argument \code{knownSegments}.
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("findLargeGaps", "default", function(chromosome=NULL, x, minLength, resolution=1L, ...) {
  # Argument 'x':
  x <- Arguments$getNumerics(x);
  nbrOfLoci <- length(x);

  # Argument 'chromosome':
  if (!is.null(chromosome)) {
    disallow <- c("Inf");
    chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), disallow=disallow, length=c(nbrOfLoci, nbrOfLoci));
  }

  # Argument 'minLength':
  minLength <- Arguments$getNumeric(minLength, range=c(0,Inf));

  # Argument 'resolution':
  resolution <- Arguments$getNumeric(resolution, range=c(0,Inf));
  if (resolution >= minLength) {
    throw(sprintf("Cannot identify large gaps. Argument 'resolution' (=%g) is not strictly smaller than 'minLength' (=%g).", resolution, minLength))
  }

  if (!is.null(chromosome)) {
    allChromosomes <- sort(unique(chromosome));
    nbrOfChromosomes <-  length(allChromosomes);

    xEmpty <- vector(mode(x), length=0L)
    gaps <- data.frame(chromosome=integer(0L), start=xEmpty, end=xEmpty);
    for (cc in seq(along=allChromosomes)) {
      chr <- allChromosomes[cc];
      idxs <- which(chromosome == chr);
      chromosomeCC <- chromosome[idxs];
      xCC <- x[idxs];
      gapsCC <- findLargeGaps(chromosome=NULL, x=xCC, minLength=minLength, ...);
      if (nrow(gapsCC) > 0) {
        gapsCC <- cbind(chromosome=chr, gapsCC);
        gaps <- rbind(gaps, gapsCC);
      }
    } # for (cc ...)
  } else {
    x <- x[is.finite(x)];
    x <- sort(x);
    dx <- diff(x);

    isGap <- (dx >= minLength);
    idxsL <- which(isGap);
##str(list(x=x, dx=dx, isGap=isGap, idxsL=idxsL))
    xL <- x[idxsL];
    xR <- x[idxsL+1L];
##str(list(x=x, dx=dx, isGap=isGap, idxsL=idxsL, xL=xL, xR=xR))
    gaps <- data.frame(start=xL+resolution, end=xR-resolution);
    gaps$length <- gaps$end - gaps$start;
  }

  ## Sanity checks
  stopifnot(is.data.frame(gaps))
  stopifnot(all(gaps$start <= gaps$end))
  stopifnot(all(gaps$length >= 0))

  gaps;
}) # findLargeGaps()


setMethodS3("findLargeGaps", "data.frame", function(chromosome, ...) {
  data <- chromosome;
  findLargeGaps(chromosome=data$chromosome, x=data$x, ...);
}) # findLargeGaps()



###############################################################################
# HISTORY:
# 2015-04-25
# o BUG FIX: findLargeGaps() could return NULL.  Now it always returns
#   a data.frame.
# 2012-02-22
# o BUG FIX: findLargeGaps() did not handle missing values for
#   argument 'chromosome'.
# 2011-11-22
# o Added findLargeGaps().
# o Created.
###############################################################################
