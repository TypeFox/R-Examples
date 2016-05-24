###########################################################################/**
# @set "class=AbstractCBS"
# @RdocDocumentation "Restructuring AbstractCBS objects"
# @alias RestructPSCBS
#
# \description{
#   This page describes available methods for restructuring an
#   @see "AbstractCBS" object.
#
#   \itemize{
#     \item @seemethod "append"
#           - Appends one @see "AbstractCBS" to another.
#
#     \item @seemethod "extractChromosomes" /
#           @seemethod "extractChromosome"
#           - Extracts an @see "AbstractCBS" with the specified chromosomes.
#
#     \item @seemethod "extractSegments" /
#           @seemethod "extractSegment"
#           - Extracts an @see "AbstractCBS" with the specified segments.
#
#     \item @seemethod "extractRegions" /
#           @seemethod "extractRegion"
#           - Extracts an @see "AbstractCBS" with the specified regions
#             each of a certain size, where a region is defined as a
#             connected set of segments.
#
#     \item @seemethod "dropRegions" /
#           @seemethod "dropRegion"
#           - Drops specified regions and returns an @see "AbstractCBS"
#             without them.
#
#     \item @seemethod "dropChangePoint" /
#           @seemethod "mergeTwoSegments"
#           - Drops a change point by merging two neighboring segments
#             and recalculates the statistics for the merged segment
#             before returning an @see "AbstractCBS".
#
#     \item @seemethod "dropChangePoints"
#           - Drops zero or more change points
#             and recalculates the segment statistics
#             before returning an @see "AbstractCBS".
#
#     \item @seemethod "mergeThreeSegments"
#           - Merges a segment with its two flanking segments
#             and recalculates the statistics for the merged segment
#             before returning an @see "AbstractCBS".
#   }
#
#   All of the above methods are implemented for @see "CBS" and
#   @see "PairedPSCBS" objects.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################


###########################################################################/**
# @set "class=AbstractCBS"
# @RdocMethod append
#
# @title "Appends one segmentation result to another"
#
# \description{
#   @get "title",
#   where both holds segmentation results \emph{of the same sample}.
# }
#
# @synopsis
#
# \arguments{
#  \item{x, other}{The two @see "AbstractCBS" objects to be combined.}
#  \item{addSplit}{If @TRUE, a "divider" is added between chromosomes.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a object of the same class as argument \code{x}.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("append", "AbstractCBS", abstract=TRUE);



setMethodS3("renameChromosomes", "AbstractCBS", function(fit, from, to, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'from' & 'to':
  from <- Arguments$getIntegers(from, disallow=c("NaN", "Inf"));
  n <- length(from);
  to <- Arguments$getIntegers(to, disallow=c("NaN", "Inf"), length=c(n,n));


  # Nothing to do?
  if (n == 0) {
    return(fit);
  }

  data <- getLocusData(fit);
  segs <- getSegments(fit, splitters=TRUE, simplify=FALSE);
  knownSegments <- fit$params$knownSegments;

  for (cc in seq(length=n)) {
    chr <- from[cc];
    chrN <- to[cc];
    data$chromosome[data$chromosome == chr] <- chrN;
    segs$chromosome[segs$chromosome == chr] <- chrN;
    knownSegments$chromosome[knownSegments$chromosome == chr] <- chrN;
  } # for (cc ...)

  fit <- setLocusData(fit, data);
  fit <- setSegments(fit, segs);
  fit$params$knownSegments <- knownSegments;

  fit;
}, protected=TRUE) # renameChromosomes()


setMethodS3("extractChromosomes", "AbstractCBS", abstract=TRUE, protected=TRUE);


setMethodS3("extractChromosome", "AbstractCBS", function(x, chromosome, ...) {
  # To please R CMD check
  this <- x;

  # Argument 'chromosome':
  chromosome <- Arguments$getInteger(chromosome, disallow=c("NaN", "Inf"));

  extractChromosomes(this, chromosomes=chromosome, ...);
}, protected=TRUE)



setMethodS3("extractSegments", "AbstractCBS", abstract=TRUE, protected=TRUE);

setMethodS3("extractSegment", "AbstractCBS", function(this, idx, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'region':
  idx <- Arguments$getIndex(idx, max=nbrOfSegments(this, splitters=TRUE));

  extractSegments(this, idxs=idx, ...);
}, private=TRUE) # extractSegment()


setMethodS3("extractRegions", "AbstractCBS", function(this, regions, H=1, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nbrOfSegments(this, splitters=TRUE);

  # Argument 'regions':
  regions <- Arguments$getIndices(regions, max=nbrOfSegments);

  # Argument 'H':
  H <- Arguments$getInteger(H, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extract regions of a certain length");

  verbose && cat(verbose, "Left-most segments of regions to be extracted:");
  verbose && str(verbose, regions);
  verbose && cat(verbose, "Number of segments in each region: ", H);


  # Identify segments to keep
  Hs <- seq(length=H);
  regions <- regions - 1L;
  regions <- as.list(regions);
  segments <- lapply(regions, FUN=function(region) region + Hs);
  segments <- unlist(segments, use.names=FALSE);
  segments <- sort(unique(segments));
  verbose && cat(verbose, "Final set of segments to be extracted:");
  verbose && str(verbose, segments);

  res <- extractSegments(this, idxs=segments, ...);

  verbose && exit(verbose);

  res;
}, protected=TRUE) # extractRegions()



setMethodS3("extractRegion", "AbstractCBS", function(this, region, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'region':
  region <- Arguments$getIndex(region, max=nbrOfSegments(this, splitters=TRUE));

  extractRegions(this, regions=region, ...);
}, private=TRUE) # extractRegion()



###########################################################################/**
# @RdocMethod mergeTwoSegments
# @aliasmethod dropChangePoint
#
# @title "Merge two neighboring segments"
#
# \description{
#   @get "title" into one segment, which is done by dropping their
#   common change point and recalculating the segment statistics.
# }
#
# @synopsis
#
# \arguments{
#  \item{left}{An @integer specifying the segments (left, left+1)
#    to be merged.}
#  \item{update}{If @TRUE, segment statistics are updated.}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @see "AbstractCBS" of the same class with one less segment.
# }
#
# @author "HB"
#
# \seealso{
#   To merge a segment and its two flanking segments, see
#   @seemethod "mergeThreeSegments".
#   To drop regions (a connected set of segments)
#   see @seemethod "dropRegions".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("mergeTwoSegments", "AbstractCBS", abstract=TRUE, protected=TRUE);


setMethodS3("dropChangePoint", "AbstractCBS", function(fit, idx, ...) {
  # Argument 'idx':
##  max <- nbrOfChangePoints(fit, splitters=TRUE, ...);
  max <- nbrOfSegments(fit, splitters=TRUE, ...) - 1L;
  idx <- Arguments$getIndex(idx, max=max);

  mergeTwoSegments(fit, left=idx, ...);
}, protected=TRUE)




###########################################################################/**
# @RdocMethod dropChangePoints
#
# @title "Drops zero or more change points"
#
# \description{
#   @get "title", which is done by dropping one change point at the
#   time using @seemethod "dropChangePoint"
#   and recalculating the segment statistics at the end.
#
#   \emph{NOTE: This method only works if there is only one chromosome.}
# }
#
# @synopsis
#
# \arguments{
#  \item{idxs}{An @integer @vector specifying the change points to be dropped.}
#  \item{update}{If @TRUE, segment statistics are updated.}
#  \item{...}{Other arguments passed to @seemethod "dropChangePoint"
#             and @seemethod "updateMeans".}
# }
#
# \value{
#   Returns an @see "AbstractCBS" of the same class with
#   \code{length(idxs)} segments.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("dropChangePoints", "AbstractCBS", function(fit, idxs, update=TRUE, ...) {
  # Assert that there is only one chromosome
  chrs <- getChromosomes(fit);
  if (length(chrs) > 1) {
    throw("dropChangePoints() only support single-chromosome data: ", hpaste(chrs));
  }

  # Argument 'idxs':
##  max <- nbrOfChangePoints(fit, splitters=TRUE, ...);
  max <- nbrOfSegments(fit, splitters=TRUE, ...) - 1L;
  idxs <- Arguments$getIndices(idxs, max=max);


  # Drop change points one by one
  idxs <- unique(idxs);
  idxs <- sort(idxs, decreasing=TRUE);
  for (ii in seq_along(idxs)) {
    idx <- idxs[ii];
    updateI <- update && (ii == length(idxs));
    fit <- dropChangePoint(fit, idx=idx, update=updateI, ...);
  }

  # Update segment statistics?
  if (update) {
    fit <- updateMeans(fit, ...);
  }

  fit;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod mergeThreeSegments
#
# @title "Merge a segment and its two flanking segments"
#
# \description{
#   @get "title" into one segment, and recalculating the segment statistics.
# }
#
# @synopsis
#
# \arguments{
#  \item{middle}{An @integer specifying the three segments
#    (middle-1, middle, middle+1) to be merged.}
#  \item{...}{Additional arguments passed to @seemethod "mergeTwoSegments".}
# }
#
# \value{
#   Returns an @see "AbstractCBS" of the same class with two less segment.
# }
#
# @author "HB"
#
# \seealso{
#   Internally @seemethod "mergeTwoSegments" is used.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("mergeThreeSegments", "AbstractCBS", function(fit, middle, ...) {
  # Argument 'middle':
  S <- nbrOfSegments(fit, splitters=TRUE);
  middle <- Arguments$getIndex(middle, range=c(2L, S));

  # Assert that the three segments are on the same chromosome
  idxs <- middle + c(-1L, 0L, +1L);
  fitT <- extractSegments(fit, idxs);
  chrs <- getChromosomes(fitT);
  if (length(chrs) != 1L) {
    throw("Argument 'middle' specifies a segment that is at the very end of a chromosome: ", middle);
  }
  fitT <- NULL; # Not needed anymore

  fit <- mergeTwoSegments(fit, left=middle, ...);
  fit <- mergeTwoSegments(fit, left=middle-1L, ...);
  fit;
}) # mergeThreeSegments()



###########################################################################/**
# @RdocMethod dropRegions
# @aliasmethod dropRegion
#
# @title "Drops chromosomal regions (a connected set of segments)"
#
# \description{
#   @get "title" each of a certain size (number of segments).
#   \emph{None of the statistics are recalculated}.
# }
#
# @synopsis
#
# \arguments{
#  \item{regions}{An @integer @vector of length R specifying the indices
#    of the left most segment in each of the R regions to be dropped.}
#  \item{H}{A non-negative @integer specifying the size of each region,
#    i.e. the number of segments per region.}
#  \item{...}{Additional arguments passed to @seemethod "extractRegions".}
#  \item{asMissing}{If @TRUE, dropped segments are replaced by missing values,
#    otherwise they are truly dropped.}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns an @see "AbstractCBS" object of the same class with (at most)
#   R*H segments dropped.
#   If some regions overlap (share segments), then fewer than R*H segments
#   are dropped.
# }
#
# @author "HB"
#
# \seealso{
#   Internally @seemethod "extractRegions" is used.
#   See also @seemethod "dropChangePoint" and @seemethod "mergeTwoSegments".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("dropRegions", "AbstractCBS", function(this, regions, H=1, ..., asMissing=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nbrOfSegments(this, splitters=TRUE);
  # Argument 'regions':
  regions <- Arguments$getIndices(regions, max=nbrOfSegments);

  # Argument 'H':
  H <- Arguments$getInteger(H, range=c(0,Inf));

  # Argument 'asMissing':
  asMissing <- Arguments$getLogical(asMissing);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Dropping regions of a certain length");

  verbose && cat(verbose, "Left-most segments of regions to be dropped:");
  verbose && str(verbose, regions);
  verbose && cat(verbose, "Number of segments in each region: ", H);

  # Nothing to do?
  if (H == 0) {
    verbose && cat(verbose, "Nothing to do. No segments will be dropped.");
    verbose && exit(verbose);
    return(this);
  }

  # Identify segments to drop
  Hs <- seq(length=H);
  regions <- regions - 1L;
  regions <- as.list(regions);
  regions <- lapply(regions, FUN=function(region) region + Hs);
  regions <- unlist(regions, use.names=FALSE);
  regions <- sort(unique(regions));
  verbose && cat(verbose, "Final set of segments to be dropped:");
  verbose && str(verbose, regions);

  # Identify segments to keep
  allRegions <- seq(length=nbrOfSegments);
  keepSegments <- setdiff(allRegions, regions);
  verbose && cat(verbose, "Final set of segments to be kept:");
  verbose && str(verbose, keepSegments);

  dropped <- extractRegions(this, regions=regions, ...);
  res <- this;
  if (length(regions) > 0) {
    if (asMissing) {
      segs <- getSegments(res, splitters=TRUE);
      pattern <- "(chromosome|id|start|end)$";

      # TODO/AD HOC: Should be class specific /HB 2011-10-17
      pattern <- "(chromosome|id)$";
      excl <- grep(pattern, colnames(segs), ignore.case=TRUE, invert=TRUE);
      segs[regions,excl] <- NA;
      res$output <- segs;

      # TODO/AD HOC: Should be class specific /HB 2011-10-17
      for (ff in grep("segRows", names(res), ignore.case=TRUE, value=TRUE)) {
        res[[ff]][regions,] <- NA;
      }
    } else {
      res <- extractRegions(res, regions=keepSegments, ...);
    }
  }
  res$dropped <- dropped;

  # Sanity check
  if (asMissing) {
    stopifnot(nbrOfSegments(res, splitters=TRUE) == nbrOfSegments(this, splitters=TRUE));
  } else {
    stopifnot(nbrOfSegments(res, splitters=TRUE) + length(regions) == nbrOfSegments(this, splitters=TRUE));
  }

  verbose && exit(verbose);

  res;
}, protected=TRUE)


setMethodS3("dropRegion", "AbstractCBS", function(fit, region, ...) {
  # Argument 'region':
  region <- Arguments$getIndex(region);

  dropRegions(fit, regions=region, ...);
}, protected=TRUE)


setMethodS3("shiftTCN", "AbstractCBS", abstract=TRUE, protected=TRUE);




############################################################################
# HISTORY:
# 2013-04-20 [HB]
# o CLEANUP: Removed previously deprecated methods for AbstractCBS.
# 2013-03-21 [HB]
# o SPEEDUP: Made dropChangePoints() faster by only updating the segment
#   statistics/means at the very end.
# o BUG FIX: dropChangePoint[s]() for AbstractCBS would not allow to
#   drop the change points at the very end, if segmentation where done
#   with known segments/gaps and/or empty segments.
# 2012-09-13
# o Now renameChromosomes() also adjusts 'knownSegments'.
# o Added shiftTCN().
# 2012-02-27
# o Added renameChromosomes() to AbstractCBS.
# 2012-02-25
# o Added dropChangePoints() for AbstractCBS.
# 2011-11-17
# o FIX: extractRegions() for AbstractCBS would also show verbose output.
# 2011-11-04
# o BUG FIX: extractSegment() for AbstractCBS would give an error, because
#   it called itself instead of extractSegments().
# 2011-10-21
# o Added mergeThreeSegments() to AbstractCBS.
# 2011-10-17
# o Added argument 'asMissing' to dropRegions() for AbstractCBS.
# 2011-10-14
# o Added implementation of extractRegions() for AbstractCBS, which
#   utilizes extractSegments().
# o Added abstract extractSegments() and extractSegment() for AbstractCBS.
# 2011-10-10
# o Added extractRegion()/dropRegion() and extractRegions()/dropRegions()
#   for AbstractCBS, where former ones are wrappers for the latter ones.
# o Added dropChangePoint() for AbstractCBS, which is just a
#   "name wrapper" for mergeTwoSegments().
# 2011-10-08
# o Added abstract updateMeans() for AbstractCBS.
# o Added abstract mergeTwoSegments() for AbstractCBS.
# 2011-10-02
# o Created.
############################################################################
