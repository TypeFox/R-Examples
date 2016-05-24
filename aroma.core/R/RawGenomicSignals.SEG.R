setMethodS3("extractDataForSegmentation", "RawGenomicSignals", function(this, order=TRUE, useWeights=TRUE, dropNonFinite=TRUE, dropZeroWeights=TRUE, dropWeightsIfAllEqual=TRUE, defaultChromosome=0L, defaultSampleName="Unnamed sample", ..., verbose=FALSE) {
  # This is a single-chromosome method. Assert that is the case.
  assertOneChromosome(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extracting data used by segmentation algorithms");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop loci with unknown locations?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Dropping loci with unknown locations");
  x <- this$x;
  keep <- which(is.finite(x));
  nbrOfDropped <- length(x)-length(keep);
  verbose && cat(verbose, "Number of dropped loci: ", nbrOfDropped);
  if (nbrOfDropped > 0) {
    this <- this[keep,,drop=FALSE];
  }
  # Not needed anymore
  x <- keep <- NULL;
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Order along genome?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (order) {
    verbose && enter(verbose, "Ordering along genome");
    this <- sort(this);
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  sampleName <- getBasicField(this, "fullname");
  if (is.null(sampleName)) {
    sampleName <- getFullName(this);
  }
  if (is.null(sampleName)) {
    sampleName <- defaultSampleName;
  }

  chromosome <- as.integer(this$chromosome);
  if (all(is.na(chromosome))) {
    chromosome <- defaultChromosome;
  }
  nbrOfLoci <- nbrOfLoci(this);
  verbose && cat(verbose, "Sample name: ", sampleName);
  verbose && cat(verbose, "Chromosomes: ", hpaste(sort(unique(chromosome))));
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);

  # Extracting data of interest
  data <- as.data.frame(this, translate=FALSE);
  if (!is.element("chromosome", colnames(data))) {
    data <- cbind(chromosome=chromosome, data);
  }
#  verbose && str(verbose, data);

  # Use weights, if they exists?
  hasWeights <- useWeights && (length(data$w) > 0);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop non-finite signals?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (dropNonFinite) {
    verbose && enter(verbose, "Dropping loci with non-finite signals");
    y <- data$y;
    verbose && cat(verbose, "Signals:");
    verbose && str(verbose, y);
    # Sanity check
    stopifnot(is.numeric(y));
    keep <- which(is.finite(y));
    nbrOfDropped <- nbrOfLoci-length(keep);
    verbose && cat(verbose, "Number of dropped loci: ", nbrOfDropped);
    if (nbrOfDropped > 0) {
      data <- data[keep,,drop=FALSE];
      nbrOfLoci <- nrow(data);
#      verbose && str(verbose, data);
    }
    # Not needed anymore
    keep <- NULL;
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (hasWeights && dropZeroWeights) {
    # Dropping loci with non-positive weights
    verbose && enter(verbose, "Dropping loci with non-positive weights");
    keep <- which(data$w > 0);
    nbrOfDropped <- nbrOfLoci-length(keep);
    verbose && cat(verbose, "Number of loci dropped: ", nbrOfDropped);
    if (nbrOfDropped > 0) {
      data <- data[keep,,drop=FALSE];
      nbrOfLoci <- nrow(data);
#      verbose && str(verbose, data);
    }
    # Not needed anymore
    keep <- NULL;
    verbose && exit(verbose);
  }


  if (hasWeights && dropWeightsIfAllEqual) {
    # Are all weights equal?
    verbose && enter(verbose, "Checking if all (remaining) weights are identical");
    t <- data$w - data$w[1];
    if (all(isZero(t))) {
      verbose && cat(verbose, "Dropping weights, because all weights are equal: ", data$w[1]);
      hasWeights <- FALSE;
      data$w <- NULL;
    }
    # Not needed anymore
    t <- NULL;
    verbose && exit(verbose);
  }

  attr(data, "sampleName") <- sampleName;
  verbose && str(verbose, data);

  verbose && exit(verbose);

  data;
}, protected=TRUE)

############################################################################
# HISTORY:
# 2009-05-10
# o This method supports all the segmentByNnn() methods.
# o Created.
############################################################################
