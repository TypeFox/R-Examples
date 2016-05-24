###########################################################################/**
# @RdocGeneric callSegmentationOutliers
# @alias callSegmentationOutliers.default
# @alias callSegmentationOutliers.data.frame
# @alias dropSegmentationOutliers
# @alias dropSegmentationOutliers.default
# @alias dropSegmentationOutliers.data.frame
#
# @title "Calls/drops single-locus outliers along the genome"
#
# \description{
#  @get "title" that have a signal that differ significantly from the
#  neighboring loci.
# }
#
# \usage{
#  @usage callSegmentationOutliers,default
#  @usage callSegmentationOutliers,data.frame
#  @usage dropSegmentationOutliers,default
#  @usage dropSegmentationOutliers,data.frame
# }
#
# \arguments{
#   \item{y}{A @numeric @vector of J genomic signals to be segmented.}
#   \item{chromosome}{(Optional) An @integer scalar
#       (or a @vector of length J contain a unique value).
#       Only used for annotation purposes.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{method}{A @character string specifying the method
#        used for calling outliers.}
#   \item{...}{Additional arguments passed to internal outlier
#        detection method.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   \code{callSegmentationOutliers()} returns a @logical @vector of length J.
#   \code{dropSegmentationOutliers()} returns an object of the same type
#   as argument \code{y}, where the signals for which outliers were called
#   have been set to @NA.
# }
#
# \section{Missing and non-finite values}{
#   Signals as well as genomic positions may contain missing
#   values, i.e. @NAs or @NaNs.  By definition, these cannot
#   be outliers.
# }
#
# @author "HB"
#
# \seealso{
#   Internally @see "DNAcopy::smooth.CNA" is utilized to identify
#   the outliers.
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("callSegmentationOutliers", "default", function(y, chromosome=0, x=NULL, method="DNAcopy::smooth.CNA", ..., verbose=FALSE) {
  # Local copies of DNAcopy functions
  CNA <- .use("CNA", package="DNAcopy");
  smooth.CNA <- .use("smooth.CNA", package="DNAcopy");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  disallow <- c("Inf");
  y <- Arguments$getDoubles(y, disallow=disallow);
  nbrOfLoci <- length(y);

  length2 <- rep(nbrOfLoci, times=2L);

  # Argument 'chromosome':
  disallow <- c("NaN", "Inf");
  chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), disallow=disallow);
  if (length(chromosome) == 1L) {
    chromosome <- rep(chromosome, times=nbrOfLoci);
  } else {
    chromosome <- Arguments$getVector(chromosome, length=length2);
  }

  # Argument 'x':
  if (!is.null(x)) {
    disallow <- c("Inf");
    x <- Arguments$getDoubles(x, length=length2, disallow=disallow);
  }

  # Argument 'method':
  method <- match.arg(method);

  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }




  verbose && enter(verbose, "Identifying outliers");
  uChromosomes <- sort(unique(chromosome));
  nbrOfChromosomes <- length(uChromosomes);
  verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes);
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);
  verbose && cat(verbose, "Detection method: ", method);

  # Allocate result vector
  isOutlier <- logical(nbrOfLoci);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Filter missing data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying loci with non-missing data");
  keep <- (!is.na(x) & !is.na(y));
  if (!is.null(chromosome)) {
    keep <- (keep & !is.na(chromosome));
  }
  keep <- which(keep);
  chromosome <- chromosome[keep];
  x <- x[keep];
  y <- y[keep];
  nbrOfLoci <- length(x);
  verbose && cat(verbose, "Number of loci with non-missing data: ", nbrOfLoci);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For each chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  isOutlierT <- logical(nbrOfLoci);
  for (kk in seq(along=uChromosomes)) {
    chr <- uChromosomes[kk];
    verbose && enter(verbose, sprintf("Chromosome #%d ('Chr%02d') of %d", kk, chr, length(uChromosomes)));
    keepKK <- which(chromosome == chr);
    nbrOfLociKK <- length(keepKK);
    verbose && cat(verbose, "Number of loci on chromosome: ", nbrOfLociKK);

    # Extract data
    yKK <- y[keepKK];
    xKK <- x[keepKK];
    chromosomeKK <- chromosome[keepKK];

    # Order loci along chromosome
    o <- order(xKK);
    xKK <- xKK[o];
    yKK <- yKK[o];
    chromosomeKK <- chromosomeKK[o];
    keepKK <- keepKK[o];
    o <- NULL; # Not needed anymore

    # Supress all warnings, in order to avoid warnings by DNAcopy::CNA()
    # on "array has repeated maploc positions".  Ideally we should filter
    # just those out. /HB 2013-10-22
    suppressWarnings({
      dataKK <- CNA(genomdat=yKK, chrom=chromosomeKK, maploc=xKK, sampleid="y", presorted=TRUE);
    });
    chromosomeKK <- xKK <- NULL; # Not needed anymore

    yKKs <- smooth.CNA(dataKK, ...)$y;
    dataKK <- NULL; # Not needed anymore

    # Sanity check
    stopifnot(length(yKKs) == nbrOfLociKK);
    outliersKK <- which(yKKs != yKK);
    yKKs <- yKK <- NULL; # Not needed anymore

    nbrOfOutliers <- length(outliersKK);
    verbose && cat(verbose, "Number of outliers: ", nbrOfOutliers);

    outliers <- keepKK[outliersKK];
    keepKK <- outliersKK <- NULL; # Not needed anymore

    isOutlierT[outliers] <- TRUE;
    outliers <- NULL; # Not needed anymore

    verbose && exit(verbose);
  } # for (kk ...)
  chromosome <- x <- y <- NULL; # Not needed anymore

  isOutlier[keep] <- isOutlierT;
  isOutlierT <- keep <- NULL; # Not needed anymore

  nbrOfOutliers <- sum(isOutlier, na.rm=TRUE);
  verbose && cat(verbose, "Total number of outliers: ", nbrOfOutliers);

  verbose && exit(verbose);

  isOutlier;
}) # callSegmentationOutliers()


setMethodS3("callSegmentationOutliers", "data.frame", function(y, ...) {
  data <- y;

  # Get either CBS or PSCBS total CN signals.
  y <- data$y;
  if (is.null(y)) {
    y <- data$CT;
  }

  callSegmentationOutliers(y=y, chromosome=data$chromosome, x=data$x, ...);
}) # callSegmentationOutliers()


setMethodS3("dropSegmentationOutliers", "default", function(y, ...) {
  isOutlier <- callSegmentationOutliers(y, ...);
  y[isOutlier] <- NA_real_;
  isOutlier <- NULL; # Not needed anymore
  y;
})


setMethodS3("dropSegmentationOutliers", "data.frame", function(y, ...) {
  data <- y;

  isOutlier <- callSegmentationOutliers(data, ...);

  # Update either CBS or PSCBS total CN signals.
  key <- "CT";
  if (!is.element(key, colnames(data))) {
    key <- "y";
  }

  data[[key]][isOutlier] <- NA_real_;

  isOutlier <- NULL; # Not needed anymore

  data;
})


############################################################################
# HISTORY:
# 2014-02-04
# o Now retrieving local copies on DNAcopy functions up front.
# 2013-12-04
# o DOCUMENTATION: Now {call|drop}SegmentationOutliers() are documented
#   as generic functions.
# o Now {call|drop}SegmentationOutliers() drops allocated memory faster.
# o Added Rdoc for dropSegmentationOutliers().
# 2011-11-23
# o Added callSegmentationOutliers() and dropSegmentationOutliers()
#   for data frames.
# 2011-05-31
# o Now explicitly using DNAcopy::nnn() to call DNAcopy functions.
# 2010-11-27
# o Added dropSegmentationOutliers() which sets outliers to missing values.
# o Added callSegmentationOutliers(), which utilizes the detection method
#   of DNAcopy::smooth.CNA() as suggested by ABO.
# o Created.
############################################################################
