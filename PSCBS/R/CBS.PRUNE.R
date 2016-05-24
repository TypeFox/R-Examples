###########################################################################/**
# @set "class=CBS"
# @RdocMethod pruneBySdUndo
#
# @title "Prune the CBS profile by dropping change points that are too small"
#
# \description{
#  @get "title", where "too small" means that the amplitude of the
#  change points is less than a multiple of the overall standard deviation
#  of the copy-number signals.
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A @see "CBS" object.}
#   \item{rho}{A positive @double scalar specifying the number of standard
#     deviations (\code{rho*sigma}) required in order to keep a change point.
#     More change points are dropped the greater this value is.}
#   \item{sigma}{The whole-genome standard deviation of the locus-level
#     copy number signals.  The default is to calculate it from the data
#     and as done in the \pkg{DNAcopy} package.}
#   \item{...}{(Optional) Additional arguments passed to the standard
#     deviation estimator function.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "CBS" object (of the same class as \code{fit}).
# }
#
# \details{
#  This method corresponds to using the \code{undo} argument when calling
#  @see "segmentByCBS", which in turn corresponds to using the
#  \code{undo.splits="sdundo"} and \code{undo.SD} of the underlying
#  @see "DNAcopy::segment" method.
# }
#
# @examples "../incl/segmentByCBS,pruneBySdUndo.Rex"
#
# @author "HB, PN"
#
# @keyword internal
#*/###########################################################################
setMethodS3("pruneBySdUndo", "CBS", function(fit, rho=3, sigma="DNAcopy", ..., verbose=FALSE) {
  # Local copies of DNAcopy functions
  DNAcopy_changepoints.sdundo <- .use("changepoints.sdundo", package="DNAcopy");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'rho':
  rho <- Arguments$getDouble(rho, range=c(0,Inf));

  # Argument 'sigma':
  if (is.character(sigma)) {
    sigma <- estimateStandardDeviation(fit, method=sigma, ...);
  }
  sigma <- Arguments$getDouble(sigma, range=c(0,Inf), disallow=c("NA", "NaN", "Inf"));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Pruning segments by standard deviation");

  # Check if locus weights are available
  data <- getLocusData(fit);
  hasWeights <- !is.null(data$w);
  # Not needed anymore
  data <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Prune chromosome by chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chromosomes <- getChromosomes(fit);
  nbrOfChromosomes <- length(chromosomes);

  fitList <- vector("list", length=nbrOfChromosomes);
  for (cc in seq(length=nbrOfChromosomes)) {
    chr <- chromosomes[cc];
    verbose && enter(verbose, sprintf("Chromosome #%d ('Chr%s') of %d",
                                            cc, chr, length(chromosomes)));

    # Extract this chromosome
    fitT <- extractChromosome(fit, chromosome=chr);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get segmentation data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    data <- getLocusData(fitT);
    segs <- getSegments(fitT);
    segRows <- fitT$segRows;
    nbrOfSegs <- nrow(segRows);
    verbose && cat(verbose, "Number of segments (before): ", nbrOfSegs);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Drop missing values
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    y <- data$y;

    # Label data points by their segment index
    segId <- rep(NA_integer_, times=max(segRows[,2], na.rm=TRUE));
    for (rr in 1:nbrOfSegs) {
      segRow <- unlist(segRows[rr,], use.names=FALSE);
      idxs <- segRow[1]:segRow[2];
      segId[idxs] <- rr;
    }

    # Drop missing value
    keep <- !is.na(y);
    if (hasWeights) {
      w <- data$w;
      keep <- keep & !is.na(w);
    }
    units <- which(keep);
    y <- y[units];
    segId <- segId[units];
    if (hasWeights) {
      w <- w[units];
    }

    # Update 'segRows' accordingly
    for (rr in 1:nbrOfSegs) {
      startStop <- range(which(segId == rr));
      segRows[rr,1] <- startStop[1];
      segRows[rr,2] <- startStop[2];
    }
    # Not needed anymore
    segId <- startStop <- NULL;

    # Sanity check
    stopifnot(max(segRows[,2]) <= length(y));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Prune change points
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    segLengths <- segRows[,2] - segRows[,1] + 1L;
    segLengthsP <- DNAcopy_changepoints.sdundo(genomdat=y,
                         lseg=segLengths, trimmed.SD=sigma, change.SD=rho);
    segLengthsP <- as.integer(segLengthsP);
    nbrOfSegsP <- length(segLengthsP);
    verbose && cat(verbose, "Number of segments (after): ", nbrOfSegsP);

    nbrOfPrunedSegs <- nbrOfSegs-nbrOfSegsP;
    verbose && cat(verbose, "Number of segments dropped: ", nbrOfPrunedSegs);

    # No segments pruned?
    if (nbrOfPrunedSegs == 0) {
      # Sanity check
      stopifnot(identical(segLengthsP, segLengths));

      fitList[[cc]] <- fitT;
      verbose && cat(verbose, "Nothing to changed. Skipping.");
#      verbose && exit(verbose);
#      next;
    }

    # Setup new 'segRows'
    endRow <- cumsum(segLengthsP);
    n <- length(endRow);
    segRowsP <- data.frame(startRow=c(1L, endRow[-n]+1L), endRow=endRow);


    # Expand to units with also missing values
    segRowsP[,1] <- units[segRowsP[,1]];
    segRowsP[,2] <- units[segRowsP[,2]];

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Create stub for a segment table
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    idxs <- seq(length=nbrOfSegsP);
    segsP <- segs[idxs,];

    # Sanity checks
    if (nbrOfPrunedSegs == 0) {
      segRows <- fitT$segRows;
      stopifnot(all.equal(segRowsP, segRows, check.attributes=FALSE));
      stopifnot(all.equal(segsP, segs, check.attributes=FALSE));
    }

    fitT$output <- segsP;
    fitT$segRows <- segRowsP;

    fitList[[cc]] <- fitT;

    verbose && exit(verbose);
  } # for (cc ...)

  fitP <- Reduce(append, fitList);

  verbose && enter(verbose, "Updating segment means and boundaries");
  fitP <- updateBoundaries(fitP, verbose=less(verbose, 50));
  fitP <- updateMeans(fitP, verbose=less(verbose, 50));
  verbose && exit(verbose);

  nbrOfSegs <- nbrOfSegments(fit);
  nbrOfSegsP <- nbrOfSegments(fitP);
  nbrOfPrunedSegs <- nbrOfSegs-nbrOfSegsP;
  verbose && cat(verbose, "Number of segments (before): ", nbrOfSegs);
  verbose && cat(verbose, "Number of segments (after): ", nbrOfSegsP);
  verbose && cat(verbose, "Number of segments dropped: ", nbrOfPrunedSegs);

  verbose && exit(verbose);

  fitP;
}) # pruneBySdUndo()


setMethodS3("seqOfSegmentsByDP", "CBS", function(fit, by=c("y"), ...) {
  NextMethod("seqOfSegmentsByDP", by=by);
})


############################################################################
# HISTORY:
# 2012-09-13
# o Added seqOfSegmentsByDP() for CBS.
# 2011-12-06
# o BUG FIX: pruneBySdUndo() for CBS did not work with more than one array.
# 2011-11-16
# o Added Rdoc comments.
# 2011-11-15
# o Added pruneBySdUndo() for CBS.
# o Created.
############################################################################
