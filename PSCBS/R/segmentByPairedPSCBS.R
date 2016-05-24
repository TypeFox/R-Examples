###########################################################################/**
# @RdocDefault segmentByPairedPSCBS
# @alias segmentByPairedPSCBS.data.frame
# @alias segmentByPairedPSCBS.PairedPSCBS
# @alias segmentByPairedPSCBS
#
# @title "Segment total copy numbers and allele B fractions using the Paired PSCBS method"
#
# \description{
#  @get "title" [1].
#  This method requires matched normals.
#  This is a low-level segmentation method.
#  It is intended to be applied to one tumor-normal sample at the time.
# }
#
# @synopsis
#
# \arguments{
#   \item{CT}{A @numeric @vector of J tumor total copy number (TCN)
#        ratios in [0,+@Inf) (due to noise, small negative values are
#        also allowed).  The TCN ratios are typically scaled such that
#        copy-neutral diploid loci have a mean of two.}
#   \item{thetaT, thetaN}{(alternative) As an alternative to specifying
#        tumor TCN \emph{ratios} relative to the match normal by
#        argument \code{CT}, on may specify total tumor and normal
#        signals seperately, in which case the TCN ratios \code{CT} are
#        calculated as \eqn{CT = 2*thetaT/thetaN}.}
#   \item{betaT}{A @numeric @vector of J tumor allele B fractions (BAFs)
#        in [0,1] (due to noise, values may be slightly outside as well)
#        or @NA for non-polymorphic loci.}
#   \item{betaN}{A @numeric @vector of J matched normal BAFs in [0,1]
#        (due to noise, values may be slightly outside as well) or @NA
#        for non-polymorphic loci.}
#   \item{muN}{An optional @numeric @vector of J genotype calls in
#        \{0,1/2,1\} for AA, AB, and BB, respectively,
#        and @NA for non-polymorphic loci.
#        If not given, they are estimated from the normal BAFs using
#        @see "aroma.light::callNaiveGenotypes" as described in [2].}
#   \item{rho}{(alternative to \code{betaT} and \code{betaN}/\code{muN})
#        A @numeric @vector of J decrease-of-heterozygosity signals (DHs)
#        in [0,1] (due to noise, values may be slightly larger than one
#        as well).  By definition, DH should be @NA for homozygous loci
#        and for non-polymorphic loci.}
#   \item{chromosome}{(Optional) An @integer scalar (or a @vector of length J),
#        which can be used to specify which chromosome each locus belongs to
#        in case multiple chromosomes are segments.
#        This argument is also used for annotation purposes.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{alphaTCN, alphaDH}{The significance levels for segmenting total
#        copy numbers (TCNs) and decrease-in-heterozygosity signals (DHs),
#        respectively.}
#   \item{undoTCN, undoDH}{Non-negative @numerics.  If greater than 0,
#        then a cleanup of segmentions post segmentation is done.
#        See argument \code{undo} of @see "segmentByCBS" for more
#        details.}
#   \item{avgTCN, avgDH}{A @character string specifying how to calculating
#         segment mean levels \emph{after} change points have been
#         identified.}
#   \item{...}{Additional arguments passed to @see "segmentByCBS".}
#   \item{flavor}{A @character specifying what type of segmentation and
#     calling algorithm to be used.}
#   \item{tbn}{If @TRUE, \code{betaT} is normalized before segmentation
#     using the TumorBoost method [2], otherwise not.}
#   \item{preserveScale}{Passed to @see "aroma.light::normalizeTumorBoost",
#     which is only called if \code{tbn} is @TRUE.}
#   \item{joinSegments}{If @TRUE, there are no gaps between neighboring
#     segments.
#     If @FALSE, the boundaries of a segment are defined by the support
#     that the loci in the segments provides, i.e. there exist a locus
#     at each end point of each segment.  This also means that there
#     is a gap between any neighboring segments, unless the change point
#     is in the middle of multiple loci with the same position.
#     The latter is what \code{DNAcopy::segment()} returns.
#   }
#   \item{knownSegments}{Optional @data.frame specifying
#     \emph{non-overlapping} known segments.  These segments must
#     not share loci.  See @see "findLargeGaps" and @see "gapsToSegments".}
#   \item{dropMissingCT}{If @TRUE, loci for which 'CT' is missing
#     are dropped, otherwise not.}
#   \item{seed}{An (optional) @integer specifying the random seed to be
#     set before calling the segmentation method.  The random seed is
#     set to its original state when exiting.  If @NULL, it is not set.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the segmentation results as a @see "PairedPSCBS" object.
# }
#
# \details{
#   Internally @see "segmentByCBS" is used for segmentation.
#   The Paired PSCBS segmentation method does \emph{not} support weights.
# }
#
# \section{Reproducibility}{
#   The "DNAcopy::segment" implementation of CBS uses approximation
#   through random sampling for some estimates.  Because of this,
#   repeated calls using the same signals may result in slightly
#   different results, unless the random seed is set/fixed.
# }
#
# \section{Whole-genome segmentation is preferred}{
#   Although it is possible to segment each chromosome independently
#   using Paired PSCBS, we strongly recommend to segment whole-genome
#   (TCN,BAF) data at once.  The reason for this is that downstream
#   CN-state calling methods, such as the AB and the LOH callers,
#   performs much better on whole-genome data.  In fact, they may
#   fail to provide valid calls if done chromsome by chromosome.
# }
#
# \section{Missing and non-finite values}{
#   The total copy number signals as well as any optional positions
#   must not contain missing values, i.e. @NAs or @NaNs.
#   If there are any, an informative error is thrown.
#   Allele B fractions may contain missing values, because such are
#   interpreted as representing non-polymorphic loci.
#
#   None of the input signals may have infinite values, i.e. -@Inf or +@Inf.
#   If so, an informative error is thrown.
# }
#
# \section{Paired PSCBS with only genotypes}{
#   If allele B fractions for the matched normal (\code{betaN}) are
#   not available, but genotypes (\code{muN}) are, then it is possible
#   to run a version of Paired PSCBS where TumorBoost normalization
#   of the tumor allele B fractions is skipped.  In order for this
#   to work, argument \code{tbn} must be set to @FALSE.
# }
#
# @examples "../incl/segmentByPairedPSCBS.Rex"
#
# @author "HB"
#
# \references{
#  [1] @include "../incl/OlshenA_etal_2011.Rd" \cr
#  [2] @include "../incl/BengtssonH_etal_2010.Rd" \cr
# }
#
# \seealso{
#   Internally, @see "aroma.light::callNaiveGenotypes" is used to
#   call naive genotypes, @see "aroma.light::normalizeTumorBoost" is
#   used for TumorBoost normalization, and @see "segmentByCBS" is used
#   to segment TCN and DH separately.
#
#   To segment tumor total copy numbers and allele B fractions
#   \emph{without} a matched normal, see @see "segmentByNonPairedPSCBS".
#
#   To segment total copy-numbers, or any other unimodal signals,
#   see @see "segmentByCBS".
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("segmentByPairedPSCBS", "default", function(CT, thetaT=NULL, thetaN=NULL, betaT=NULL, betaN=NULL, muN=NULL, rho=NULL, chromosome=0, x=NULL, alphaTCN=0.009, alphaDH=0.001, undoTCN=0, undoDH=0, ..., avgTCN=c("mean", "median"), avgDH=c("mean", "median"), flavor=c("tcn&dh", "tcn,dh", "sqrt(tcn),dh", "sqrt(tcn)&dh", "tcn"), tbn=is.null(rho), preserveScale=getOption("PSCBS/preserveScale", FALSE), joinSegments=TRUE, knownSegments=NULL, dropMissingCT=TRUE, seed=NULL, verbose=FALSE) {
  # WORKAROUND: If Hmisc is loaded after R.utils, it provides a buggy
  # capitalize() that overrides the one we want to use. Until PSCBS
  # gets a namespace, we do the following workaround. /HB 2011-07-14
  capitalize <- R.utils::capitalize;

  # To please R CMD check
  index <- NULL; rm(list="index");

  # Settings for sanity checks
  tol <- getOption("PSCBS/sanityChecks/tolerance", 0.0005);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'thetaT' & 'thetaN':
  if (!is.null(thetaT) && !is.null(thetaN)) {
    thetaT <- Arguments$getDoubles(thetaT, disallow=disallow);
    nbrOfLoci <- length(thetaT);
    length2 <- rep(nbrOfLoci, times=2L);
    thetaN <- Arguments$getDoubles(thetaN, length=length2, disallow=disallow);
    CT <- 2 * thetaT / thetaN;
  } else if (!is.null(thetaT) || !is.null(thetaN)) {
    throw("Either argument 'CT' needs to be specified or *both* of arguments 'thetaT' and 'thetaN'");
  }

  # Argument 'CT':
  disallow <- c("Inf");
  CT <- Arguments$getDoubles(CT, disallow=disallow);
  nbrOfLoci <- length(CT);
  length2 <- rep(nbrOfLoci, times=2L);


  # Argument 'betaT':
  if (!is.null(betaT)) {
    betaT <- Arguments$getDoubles(betaT, length=length2, disallow="Inf")
  }

  # Argument 'betaN':
  if (!is.null(betaN)) {
    betaN <- Arguments$getDoubles(betaN, length=length2, disallow="Inf");
  }

  # Argument 'muN':
  if (!is.null(muN)) {
    muN <- Arguments$getDoubles(muN, length=length2, range=c(0,1), disallow="Inf");
    if (all(is.na(muN)) == nbrOfLoci) {
      throw(sprintf("All genotypes ('muN') are NAs: %d (100%%) out of %d", nbrOfLoci, nbrOfLoci));
    }
  }

  # Argument 'rho':
  if (!is.null(rho)) {
    rho <- Arguments$getDoubles(rho, range=c(0,Inf), length=length2, disallow="Inf")
  }

  if (is.null(muN)) {
    if (is.null(betaN) && is.null(rho)) {
      throw("If argument 'muN' is not given, then either 'betaN' or 'rho' must be.")
    }
  }

  # Argument 'tbn':
  tbn <- Arguments$getLogical(tbn);
  if (!is.null(tbn)) {
    if (tbn) {
      if (is.null(betaT)) {
        throw("Cannot do TumorBoost normalization (tbn=TRUE) without tumor BAFs ('betaT').")
      }
      if (is.null(betaN)) {
        throw("Cannot do TumorBoost normalization (tbn=TRUE) with normal BAFs ('betaN').")
      }
    }
  }

  # Argument 'preserveScale':
  if (tbn && missing(preserveScale)) {
    if (!is.element("PSCBS/preserveScale", names(options()))) {
      warning("Argument 'preserveScale' for segmentByPairedPSCBS() now defaults to FALSE. Prior to PSCBS v0.50.0 (October 2015) the default was TRUE.  To avoid this warning, explicitly specify this argument when calling segmentByPairedPSCBS() or make sure to set option 'PSCBS/preserveScale' to either TRUE or FALSE.  This warning will be removed in a future version.");
    }
  }
  preserveScale <- Arguments$getLogical(preserveScale);

  # Argument 'chromosome':
  if (is.null(chromosome)) {
    chromosome <- 0L;
  } else {
    disallow <- c("Inf");
    chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), disallow=disallow);
    if (length(chromosome) > 1) {
      chromosome <- Arguments$getIntegers(chromosome, length=length2, disallow=disallow);
    }
  }

  # Argument 'x':
  if (is.null(x)) {
    x <- seq(length=nbrOfLoci);
  } else {
    disallow <- c("Inf");
    x <- Arguments$getDoubles(x, length=length2, disallow=disallow);
  }

  # Argument 'alphaTCN':
  alphaTCN <- Arguments$getDouble(alphaTCN, range=c(0,1));

  # Argument 'alphaDH':
  alphaDH <- Arguments$getDouble(alphaDH, range=c(0,1));

  # Argument 'undoTCN':
  undoTCN <- Arguments$getDouble(undoTCN, range=c(0,Inf));

  # Argument 'undoDH':
  undoDH <- Arguments$getDouble(undoDH, range=c(0,Inf));

  # Argument 'avgTCN' & 'avgDH':
  avgTCN <- match.arg(avgTCN);
  avgDH <- match.arg(avgDH);

  # Argument 'flavor':
  flavor <- match.arg(flavor);
  knownFlavors <- eval(formals(segmentByPairedPSCBS.default)$flavor);
  if (!is.element(flavor, knownFlavors)) {
    throw("Segmentation flavor is not among the supported ones (", paste(sprintf("\"%s\"", knownFlavors), collapse=", "), "): ", flavor);
  }

  # Argument 'joinSegments':
  joinSegments <- Arguments$getLogical(joinSegments);

  # Argument 'knownSegments':
  if (is.null(knownSegments)) {
    knownSegments <- data.frame(chromosome=integer(0), start=integer(0), end=integer(0));
  } else {
    if (!joinSegments) {
##      warning("Argument 'knownSegments' should only be specified if argument 'joinSegments' is TRUE.");
    }
  }

  if (!is.data.frame(knownSegments)) {
    throw("Argument 'knownSegments' is not a data.frame: ", class(knownSegments)[1]);
  }

  if (!all(is.element(c("chromosome", "start", "end"), colnames(knownSegments)))) {
    throw("Argument 'knownSegments' does not have the required column names: ", hpaste(colnames(knownSegments)));
  }

  # Argument 'dropMissingCT':
  dropMissingCT <- Arguments$getLogical(dropMissingCT);
  if (!dropMissingCT) {
    if (is.element(flavor, c("tcn&dh", "sqrt(tcn)&dh"))) {
      throw("Missing values in 'CT' are (currently) not supported by the chosen 'flavor': ", flavor);
    }
  }


  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getIntegers(seed);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Segmenting paired tumor-normal signals using Paired PSCBS");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call genotypes?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Are genotype calls muN missing and can they be called?
  if (is.null(muN) && !is.null(betaN)) {
    verbose && enter(verbose, "Calling genotypes from normal allele B fractions");
    verbose && str(verbose, betaN);
    callNaiveGenotypes <- .use("callNaiveGenotypes", package="aroma.light");
    muN <- callNaiveGenotypes(betaN, censorAt=c(0,1));
    verbose && cat(verbose, "Called genotypes:");
    verbose && str(verbose, muN);
    verbose && print(verbose, table(muN));
    # Assert proper calls
    muN <- Arguments$getDoubles(muN, length=length2, range=c(0,1), disallow="Inf");
    # Sanity check
    if (all(is.na(muN))) {
      throw(sprintf("All genotypes ('muN') called from the normal allele B fractions ('betaN') are NAs: %d (100%%) out of %d", nbrOfLoci, nbrOfLoci));
    }
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize betaT using betaN (TumorBoost normalization)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (tbn) {
    verbose && enter(verbose, "Normalizing betaT using betaN (TumorBoost)");
    normalizeTumorBoost <- .use("normalizeTumorBoost", package="aroma.light");
    betaTN <- normalizeTumorBoost(betaT=betaT, betaN=betaN, muN=muN, preserveScale=preserveScale);
    verbose && cat(verbose, "Normalized BAFs:");
    verbose && str(verbose, betaTN);

    # Assert that no missing values where introduced
    keep <- (is.finite(betaT) & is.finite(betaN) & is.finite(muN));
    if (anyNA(betaTN[keep])) {
      throw("Internal error: normalizeTumorBoost() introduced missing values.");
    }
    # Not needed anymore
    keep <- NULL;
    verbose && exit(verbose);
  } else {
    betaTN <- betaT;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setup up data");
  data <- data.frame(chromosome=chromosome, x=x, CT=CT)
  if (!is.null(thetaT)) {
    data$thetaT <- thetaT
    data$thetaN <- thetaN
  }
  if (!is.null(betaT)) data$betaT <- betaT
  if (!is.null(betaTN)) data$betaTN <- betaTN
  if (!is.null(betaN)) data$betaN <- betaN
  if (!is.null(muN)) data$muN <- muN
  if (!is.null(rho)) data$rho <- rho
  verbose && str(verbose, data)
  # Not needed anymore
  chromosome <- x <- CT <- thetaT <- thetaN <- betaT <- betaTN <- betaN <- muN <- rho <- NULL

  # Sanity check
  stopifnot(nrow(data) == nbrOfLoci);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop data points without known genomic positions, because that
  # is what DNAcopy::CNA() will do otherwise.  At the end, we will
  # undo this such that the returned 'data' object is complete.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ok <- (!is.na(data$chromosome) & !is.na(data$x));
  if (any(!ok)) {
    verbose && enter(verbose, "Dropping loci with unknown locations");
    verbose && cat(verbose, "Number of loci dropped: ", sum(!ok));
    data <- data[ok,,drop=FALSE];
    nbrOfLoci <- nrow(data);
    verbose && exit(verbose);
  }
  ok <- NULL; # Not needed anymore

  # Sanity check
  stopifnot(nrow(data) == nbrOfLoci);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop loci for which CT is missing (regardless of betaT)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (dropMissingCT) {
    ok <- (!is.na(data$CT))
    if (any(!ok)) {
      verbose && enter(verbose, "Dropping loci for which TCNs are missing");
      verbose && cat(verbose, "Number of loci dropped: ", sum(!ok));
      data <- data[ok,,drop=FALSE];
      nbrOfLoci <- nrow(data);
      verbose && exit(verbose);
    }
    ok <- NULL; # Not needed anymore

    # Sanity check
    stopifnot(nrow(data) == nbrOfLoci);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reorder data points along the genome, because that is what
  # DNAcopy::segment() will return.  At the end, we will undo
  # the sort such that the returned 'data' object is always in
  # the same order and number of loci as the input data.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Ordering data along genome");
  o <- order(data$chromosome, data$x, decreasing=FALSE, na.last=TRUE);
  # Any change?
  if (any(o != seq(along=o))) {
    data <- data[o,,drop=FALSE];
  }
  o <- NULL; # Not needed anymore
  verbose && str(verbose, data);
  verbose && exit(verbose);

  # Attach 'index' (guaranteed to be ordered)
  data$index <- seq(length=nrow(data));

  # Sanity check
  stopifnot(nrow(data) == nbrOfLoci);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert no missing values in (chromosome, x, CT)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  ok <- (!is.na(data$chromosome) & !is.na(data$x));
  if (!all(ok)) {
    throw("INTERNAL ERROR: Detected (chromosome, x) with missing values also after filtering.");
  }

  # Sanity check
  if (dropMissingCT) {
    ok <- (!is.na(data$CT));
    if (!all(ok)) {
      throw("INTERNAL ERROR: Detected TCN with missing values also after filtering.");
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Multiple chromosomes?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify all chromosomes, excluding missing values
  chromosomes <- sort(unique(data$chromosome), na.last=NA);
  nbrOfChromosomes <- length(chromosomes);
  if (nbrOfChromosomes > 1) {
    verbose && enter(verbose, "Segmenting multiple chromosomes");
    verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes);

    # Generate random seeds?
    seeds <- NULL
    if (!is.null(seed)) {
      randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
      verbose && printf(verbose, "Random seed temporarily set (seed=c(%s), kind=\"L'Ecuyer-CMRG\")\n", paste(seed, collapse=", "))
      seeds <- randomSeed("advance", n=nbrOfChromosomes)
      verbose && printf(verbose, "Produced %d seeds from this stream for future usage\n", length(seeds))
      randomSeed("reset")
    }

    fitList <- listenv()
    for (kk in seq(length=nbrOfChromosomes)) {
      chromosomeKK <- chromosomes[kk];
      chrTag <- sprintf("Chr%02d", chromosomeKK);
      verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", kk, chrTag, nbrOfChromosomes));

      seedKK <- seeds[[kk]]

      # Extract subset of data and parameters for this chromosome
      dataKK <- subset(data, chromosome == chromosomeKK);
      verbose && str(verbose, dataKK);
      fields <- attachLocally(dataKK, fields=c("CT", "thetaT", "thetaN", "betaT", "betaTN", "betaN", "muN", "rho", "chromosome", "x"));
      dataKK <- NULL; # Not needed anymore

      knownSegmentsKK <- NULL;
      if (!is.null(knownSegments)) {
        knownSegmentsKK <- subset(knownSegments, chromosome == chromosomeKK);
        verbose && cat(verbose, "Known segments:");
        verbose && print(verbose, knownSegmentsKK);
      }

      fitList[[chrTag]] %<=% {
        fit <- segmentByPairedPSCBS(CT=CT, thetaT=thetaT, thetaN=thetaN,
                  betaT=betaTN, betaN=betaN, muN=muN, rho=rho,
                  chromosome=chromosome, x=x,
                  tbn=FALSE, joinSegments=joinSegments,
                  knownSegments=knownSegmentsKK,
                  alphaTCN=alphaTCN, alphaDH=alphaDH,
                  undoTCN=undoTCN, undoDH=undoDH,
                  avgTCN=avgTCN, avgDH=avgDH,
                  flavor=flavor,
                  ...,
                  seed=seedKK,
                  verbose=verbose)

        # Sanity checks
        if (nrow(knownSegmentsKK) == 0) {
          stopifnot(nrow(fit$data) == length(CT))
          stopifnot(all.equal(fit$data$CT, CT))
          stopifnot(all.equal(fit$data$muN, muN))
        }

        # Update betaT (which is otherwise equals betaTN)
        fit$data$betaT <- betaT

        verbose && print(verbose, head(as.data.frame(fit)))
        verbose && print(verbose, tail(as.data.frame(fit)))

        fit
      } ## fitList[[chrTag]] <- ...

      rm(list=fields) # Not needed anymore
      verbose && exit(verbose);
    } # for (kk ...)

    verbose && enter(verbose, "Merging (independently) segmented chromosome");
    fitList <- as.list(fitList)
    fit <- Reduce(append, fitList);
    fitList <- NULL; # Not needed anymore
    verbose && str(verbose, fit);
    verbose && exit(verbose);

    # Update parameters that otherwise may be incorrect
    fit$params$tbn <- tbn;
    fit$params$seed <- seed;

    segs <- as.data.frame(fit);
    if (nrow(segs) < 6) {
      verbose && print(verbose, segs);
    } else {
      verbose && print(verbose, head(segs));
      verbose && print(verbose, tail(segs));
    }

    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Return results
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return(fit);
  } # if (nbrOfChromosomes > 1)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset 'knownSegments'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Keeping only current chromosome for 'knownSegments'");

  currChromosome <- data$chromosome[1];
  verbose && cat(verbose, "Chromosome: ", currChromosome);

  knownSegments <- subset(knownSegments, chromosome == currChromosome);
  nbrOfSegments <- nrow(knownSegments);

  verbose && cat(verbose, "Known segments for this chromosome:");
  verbose && print(verbose, knownSegments);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity checks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Here 'knownSegments' should specify at most a single chromosome
  uChromosomes <- sort(unique(knownSegments$chromosome));
  if (length(uChromosomes) > 1) {
    throw("INTERNAL ERROR: Argument 'knownSegments' specifies more than one chromosome: ", hpaste(uChromosomes));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert no missing values in (chromosome, x, CT)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  ok <- (!is.na(data$chromosome) & !is.na(data$x));
  if (!all(ok)) {
    throw("INTERNAL ERROR: Detected (chromosome, x) with missing values also after filtering.");
  }

  # Sanity check
  if (dropMissingCT) {
    ok <- (!is.na(data$CT));
    if (!all(ok)) {
      throw("INTERNAL ERROR: Detected TCN with missing values also after filtering.");
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup input data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "alphaTCN: ", alphaTCN);
  verbose && cat(verbose, "alphaDH: ", alphaDH);
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate decrease-of-heterozygosity signals (DHs)?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(data$rho)) {
    verbose && enter(verbose, "Calculating DHs")
    # SNPs are identifies as those loci that have non-missing 'betaTN' & 'muN'
    isSnp <- (!is.na(data$betaTN) & !is.na(data$muN))
    nbrOfSnps <- sum(isSnp)
    verbose && cat(verbose, "Number of SNPs: ", nbrOfSnps)

    # DH is by definition only defined for heterozygous SNPs.
    # For simplicity, we set it to be NA for non-heterozygous loci.
    isHet <- isSnp & (data$muN == 1/2)
    verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                       sum(isHet), 100*sum(isHet)/nbrOfSnps)
    rho <- rep(NA_real_, length=nbrOfLoci)
    rho[isHet] <- 2*abs(data$betaTN[isHet]-1/2)
    verbose && cat(verbose, "Normalized DHs:")
    verbose && str(verbose, rho)
    data$rho <- rho
    isSnp <- isHet <- rho <- NULL # Not needed anymore
    verbose && exit(verbose)
  }
  ## Sanity check
  stopifnot(!is.null(data$rho))



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate random seeds?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  seeds <- NULL
  if (!is.null(seed)) {
    randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
    verbose && printf(verbose, "Random seed temporarily set (seed=c(%s), kind=\"L'Ecuyer-CMRG\")\n", paste(seed, collapse=", "))
    seeds <- randomSeed("advance", n=2L) ## For TCN and DH
    names(seeds) <- c("TCN", "DH")
    verbose && printf(verbose, "Produced %d seeds from this stream for future usage\n", length(seeds))
    randomSeed("reset")
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1a. Identification of change points in total copy numbers
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identification of change points by total copy numbers");

  fields <- attachLocally(data, fields=c("CT", "thetaT", "thetaN", "chromosome", "x", "index"));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert no missing values in (chromosome, x, CT)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  ok <- (!is.na(data$chromosome) & !is.na(data$x));
  if (!all(ok)) {
    throw("INTERNAL ERROR: Detected (chromosome, x) with missing values also after filtering.");
  }

  # Sanity check
  if (dropMissingCT) {
    ok <- (!is.na(data$CT));
    if (!all(ok)) {
      throw("INTERNAL ERROR: Detected CT with missing values also after filtering.");
    }
  }


  # Segment TCN ratios
  # Calculate tumor-normal TCN ratios?
  fit <- segmentByCBS(CT,
                      chromosome=chromosome, x=x, index=index,
                      joinSegments=joinSegments,
                      knownSegments=knownSegments,
                      alpha=alphaTCN, undo=undoTCN, ...,
                      seed=seeds[["TCN"]],
                      verbose=verbose);
  verbose && str(verbose, fit);

  rm(list=fields); # Not needed anymore

  # Sanity check
  if (nrow(knownSegments) == 0) {
    stopifnot(nrow(fit$data) == nrow(data));
    stopifnot(all(fit$data$chromosome == data$chromosome));
    stopifnot(all(fit$data$x == data$x));
    stopifnot(all(fit$data$index == data$index));
    stopifnot(all.equal(fit$data$y, data$CT));
  }

  tcnSegments <- fit$output;
  tcnSegRows <- fit$segRows;
  fit <- NULL; # Not needed anymore

  # Sanity checks
  stopifnot(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE));
  stopifnot(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE));

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1b. Restructure TCN segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Restructure TCN segmentation results");
  # Drop dummy columns
  keep <- setdiff(colnames(tcnSegments), c("sampleName"));
  tcnSegments <- tcnSegments[,keep,drop=FALSE];

  # Tag fields by TCN
  names <- names(tcnSegments);
  # Adding 'tcn' prefix to column names
  names <- sprintf("tcn%s", capitalize(names));
  names <- gsub("tcnChromosome", "chromosome", names, fixed=TRUE);
  names(tcnSegments) <- names;
  verbose && print(verbose, tcnSegments);

  nbrOfSegs <- nrow(tcnSegments);
  verbose && cat(verbose, "Number of TCN segments: ", nbrOfSegs);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2a. Identification of additional change points using DH
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For each segment independently, segment decrease of heterozygousity (DH)
  # using CBS. By definition, only heterozygous SNPs are used.

  if (flavor == "tcn") {
    verbose && enter(verbose, "TCN-only segmentation");

    tcnSegsExpanded <- tcnSegRows;
    dhSegRows <- tcnSegRows;

    # Segments
    segs <- tcnSegments;
    segs[,"tcnId"] <- seq(length=nbrOfSegs);
    segs[,"dhId"] <- rep(1L, times=nbrOfSegs);
    segs[,c("tcnNbrOfSNPs", "tcnNbrOfHets", "dhNbrOfLoci")] <- 0L;
    segs[,"dhStart"] <- segs[,"tcnStart"];
    segs[,"dhEnd"] <- segs[,"tcnEnd"];

    # For each TCN segment...
    for (kk in seq(length=nbrOfSegs)) {
      tcnId <- kk;

      xStart <- tcnSegments[kk,"tcnStart"];
      xEnd <- tcnSegments[kk,"tcnEnd"];
      regionTag <- sprintf("[%10g,%10g]", xStart, xEnd);
      verbose && enter(verbose, sprintf("Total CN segment #%d (%s) of %d", kk, regionTag, nbrOfSegs));

      # Empty segment?
      rowStart <- tcnSegRows[kk,1];
      rowEnd <- tcnSegRows[kk,2];

      # Empty segment or a segment separator?
      isEmptySegment <- (is.na(rowStart) && is.na(rowEnd));

      # Nothing to do?
      if (isEmptySegment) {
        verbose && exit(verbose);
        next;
      }

      nbrOfTCNLociKK <- tcnSegments[kk,"tcnNbrOfLoci"];
      verbose && cat(verbose, "Number of TCN loci in segment: ", nbrOfTCNLociKK);
      rows <- seq(from=rowStart, length=nbrOfTCNLociKK);
      dataKK <- data[rows,,drop=FALSE];
      nbrOfLociKK <- nrow(dataKK);

      verbose && cat(verbose, "Locus data for TCN segment:");
      verbose && str(verbose, dataKK);

      verbose && cat(verbose, "Number of loci: ", nbrOfLociKK);
      hasDH <- !is.null(dataKK$rho)
      if (hasDH) {
        isSnpKK <- !is.na(dataKK$rho)
        isHetsKK <- (isSnpKK & (dataKK$rho > 0))
      } else {
        isSnpKK <- !is.na(dataKK$muN)
        isHetsKK <- (isSnpKK & (dataKK$muN == 1/2))
      }
      nbrOfSnpsKK <- sum(isSnpKK)
      nbrOfHetsKK <- sum(isHetsKK)
      verbose && printf(verbose, "Number of SNPs: %d (%.2f%%)\n",
                                    nbrOfSnpsKK, 100*nbrOfSnpsKK/nbrOfLociKK);

      verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                    nbrOfHetsKK, 100*nbrOfHetsKK/nbrOfSnpsKK);

      segs[kk,"tcnNbrOfSNPs"] <- nbrOfSnpsKK;
      segs[kk,"tcnNbrOfHets"] <- nbrOfHetsKK;
      segs[kk,"dhNbrOfLoci"] <- nbrOfHetsKK;

      # Adjust 'dhRows[kk,]'
      rows <- rows[isHetsKK];
      rows <- range(rows, na.rm=TRUE);
      dhSegRows[kk,] <- rows;

      # Sanity check
      if (nbrOfHetsKK > 0) {
        stopifnot(all(dhSegRows[kk,1] <= dhSegRows[kk,2], na.rm=TRUE));
      }

      # Calculate dhMean
      rhoKK <- dataKK[["rho"]][isHetsKK];
      segs[kk,"dhMean"] <- mean(rhoKK, na.rm=TRUE);

      verbose && exit(verbose);
    } # for (kk ...)

    # Reorder segmentation columns
    keys <- c("tcnId", "dhId", colnames(tcnSegments));
    keys <- c(keys, setdiff(colnames(segs), keys));
    segs <- segs[,keys];

    verbose && exit(verbose);
  } else {
    dhSegRows <- NULL;
    tcnSegsExpanded <- NULL;

    # For each TCN segment...
    segs <- vector("list", length=nbrOfSegs);
    for (kk in seq(length=nbrOfSegs)) {
      tcnId <- kk;

      xStart <- tcnSegments[kk,"tcnStart"];
      xEnd <- tcnSegments[kk,"tcnEnd"];
      regionTag <- sprintf("[%10g,%10g]", xStart, xEnd);
      verbose && enter(verbose, sprintf("Total CN segment #%d (%s) of %d", kk, regionTag, nbrOfSegs));

      # Empty segment?
      rowStart <- tcnSegRows[kk,1];
      rowEnd <- tcnSegRows[kk,2];

      # Empty segment or a segment separator?
      isEmptySegment <- (is.na(rowStart) && is.na(rowEnd));
      isSplitter <- (isEmptySegment && is.na(xStart) && is.na(xEnd));
      isEmptySegment <- (isEmptySegment & !isSplitter);

      if (isSplitter) {
        verbose && cat(verbose, "No signals to segment. Just a \"splitter\" segment. Skipping.");

        # Sanity check
        stopifnot(kk >= 1);

        # Add a splitter segment
        segT <- segs[[kk-1]];
        segT <- segT[NA_integer_,];
        keys <- colnames(tcnSegments);
        segT[,keys] <- tcnSegments[kk,keys];
        segT[,"tcnId"] <- tcnId;
        segT[,"dhId"] <- 1L;
        segT[,c("tcnNbrOfSNPs", "tcnNbrOfHets", "dhNbrOfLoci")] <- 0L;
        segT[,"dhStart"] <- xStart;
        segT[,"dhEnd"] <- xEnd;
        segs[[kk]] <- segT;
        verbose && print(verbose, segT);

        # Add a splitter to TCN and DH segment row matrix
        segRowsT <- dhSegRows[NA_integer_,];
        dhSegRows <- rbind(dhSegRows, segRowsT);

        segRowsT <- tcnSegsExpanded[NA_integer_,];
        tcnSegsExpanded <- rbind(tcnSegsExpanded, segRowsT);

        verbose && exit(verbose);
        next;
      } # if (isSplitter)


      nbrOfTCNLociKK <- tcnSegments[kk,"tcnNbrOfLoci"];
      verbose && cat(verbose, "Number of TCN loci in segment: ", nbrOfTCNLociKK);

      # Sanity check
      stopifnot(!isEmptySegment || (isEmptySegment && (nbrOfTCNLociKK == 0)));

      if (nbrOfTCNLociKK > 0) {
        # Extract locus data for TCN segment
        rows <- rowStart:rowEnd;
    ##    if (nrow(knownSegments) == 0) {
    ##      gammaT <- tcnSegments[kk,"tcnMean"];
    ##      verbose && print(verbose, all.equal(mean(dataKK$CT, na.rm=TRUE), gammaT, tolerance=tol));
    ##      stopifnot(all.equal(mean(dataKK$CT, na.rm=TRUE), gammaT, tolerance=tol));
    ##    }
      } else {
        rows <- integer(0);
      } # if (nbrOfTCNLociKK > 0)

      dataKK <- data[rows,,drop=FALSE];
      nbrOfLociKK <- nrow(dataKK);

      # Sanity check
      stopifnot(length(dataKK$CT) == nbrOfTCNLociKK);
    ##  stopifnot(sum(!is.na(dataKK$CT)) == nbrOfTCNLociKK);

      verbose && cat(verbose, "Locus data for TCN segment:");
      verbose && str(verbose, dataKK);

      verbose && cat(verbose, "Number of loci: ", nbrOfLociKK);

      hasDH <- !is.null(dataKK$rho)
      if (hasDH) {
        isSnpKK <- !is.na(dataKK$rho)
        isHetsKK <- (isSnpKK & (dataKK$rho > 0))
      } else {
        isSnpKK <- !is.na(dataKK$muN)
        isHetsKK <- (isSnpKK & (dataKK$muN == 1/2))
      }
      nbrOfSnpsKK <- sum(isSnpKK)
      nbrOfHetsKK <- sum(isHetsKK)
      verbose && printf(verbose, "Number of SNPs: %d (%.2f%%)\n",
                                    nbrOfSnpsKK, 100*nbrOfSnpsKK/nbrOfLociKK);
      verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                    nbrOfHetsKK, 100*nbrOfHetsKK/nbrOfSnpsKK);

      # Since segments in 'knownSegments' has already been used in the TCN
      # segmentation, they are not needed in the DH segmentation.
      currChromosome <- data$chromosome[1];
      verbose && cat(verbose, "Chromosome: ", currChromosome);
      knownSegmentsT <- data.frame(chromosome=currChromosome, start=xStart, end=xEnd);

      verbose && enter(verbose, "Segmenting DH signals");
      fields <- attachLocally(dataKK, fields=c("chromosome", "x", "rho", "index"));

      fit <- segmentByCBS(rho,
                          chromosome=chromosome, x=x,
                          joinSegments=joinSegments,
                          knownSegments=knownSegmentsT,
                          alpha=alphaDH, undo=undoDH, ...,
                          seed=seeds[["DH"]],
                          verbose=verbose);
      verbose && str(verbose, fit);
      dhSegments <- fit$output;
      dhSegRowsKK <- fit$segRows;

      verbose && cat(verbose, "DH segmentation (locally-indexed) rows:");
      verbose && print(verbose, dhSegRowsKK);
      verbose && str(verbose, index);

      # Remap to genome-wide indices
      for (cc in 1:2) {
        dhSegRowsKK[,cc] <- index[dhSegRowsKK[,cc]];
      }

      verbose && cat(verbose, "DH segmentation rows:");
      verbose && print(verbose, dhSegRowsKK);

      # Not needed anymore
      rm(list=fields);
      fit <- NULL;
      verbose && exit(verbose);

      # Drop dummy columns
      keep <- setdiff(colnames(dhSegments), c("sampleName", "chromosome"));
      dhSegments <- dhSegments[,keep,drop=FALSE];

      # Tag fields by DH
      names <- names(dhSegments);
      # Adding 'dh' prefix to column names
      names <- sprintf("dh%s", capitalize(names));
      names(dhSegments) <- names;

      # Special case: If there where not enough data to segment DH...
      if (nrow(dhSegments) == 0) {
        dhSegments <- dhSegments[NA_integer_,,drop=FALSE];
        dhSegRowsKK <- dhSegRowsKK[NA_integer_,,drop=FALSE];
      }

      verbose && cat(verbose, "DH segmentation table:");
      verbose && print(verbose, dhSegments);
      verbose && print(verbose, dhSegRowsKK);

      # Expand the TCN segmentation result data frame
      rows <- rep(kk, times=nrow(dhSegments));
      verbose && cat(verbose, "Rows:");
      verbose && print(verbose, rows);
      tcnSegmentsKK <- tcnSegments[rows,,drop=FALSE];
      tcnSegRowsKK <- tcnSegRows[rows,,drop=FALSE];
      # Sanity check
      stopifnot(nrow(tcnSegmentsKK) == nrow(dhSegments));
      stopifnot(nrow(tcnSegRowsKK) == nrow(dhSegments));
      stopifnot(is.na(tcnSegRowsKK[,1]) || is.na(dhSegRowsKK[,1]) || (tcnSegRowsKK[,1] <= dhSegRowsKK[,1]));
      stopifnot(is.na(tcnSegRowsKK[,2]) || is.na(dhSegRowsKK[,2]) || (dhSegRowsKK[,2] <= tcnSegRowsKK[,2]));
      verbose && cat(verbose, "TCN segmentation rows:");
      verbose && print(verbose, tcnSegRowsKK);
      stopifnot(all(tcnSegRowsKK[,1] == tcnSegRowsKK[1,1], na.rm=TRUE));
      stopifnot(all(tcnSegRowsKK[,2] == tcnSegRowsKK[1,2], na.rm=TRUE));

      verbose && cat(verbose, "TCN and DH segmentation rows:");
      verbose && print(verbose, tcnSegRowsKK);
      verbose && print(verbose, dhSegRowsKK);
      verbose && print(verbose, tcnSegsExpanded);

      # Append
      tcnSegsExpanded <- rbind(tcnSegsExpanded, tcnSegRowsKK);
      verbose && cat(verbose, "TCN segmentation (expanded) rows:");
      verbose && print(verbose, tcnSegsExpanded);
      rownames(tcnSegsExpanded) <- NULL;

      dhSegRows <- rbind(dhSegRows, dhSegRowsKK);
      rownames(dhSegRows) <- NULL;

      verbose && cat(verbose, "TCN and DH segmentation rows:");
      verbose && print(verbose, tcnSegRows);
      verbose && print(verbose, dhSegRows);
      verbose && print(verbose, tcnSegsExpanded);

      # Sanity checks
      stopifnot(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE));
      stopifnot(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE));
      stopifnot(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE));
      stopifnot(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE));
      stopifnot(all(tcnSegsExpanded[,1] <= tcnSegsExpanded[,2], na.rm=TRUE));
      stopifnot(all(tcnSegsExpanded[,1] <= dhSegRows[,1], na.rm=TRUE));
      stopifnot(all(tcnSegsExpanded[,2] >= dhSegRows[,2], na.rm=TRUE));
  ##    if (!all(tcnSegsExpanded[-nrow(tcnSegsExpanded),2] < tcnSegsExpanded[-1,1], na.rm=TRUE)) {
  ##      stopifnot(all(tcnSegsExpanded[-nrow(tcnSegsExpanded),2] < tcnSegsExpanded[-1,1], na.rm=TRUE));
  ##    }


      # Sanity check
      stopifnot(nrow(dhSegRows) == nrow(tcnSegsExpanded));

      # Append information on number of SNPs and hets in CN region
      tcnSegmentsKK <- cbind(
        tcnSegmentsKK,
        tcnNbrOfSNPs=nbrOfSnpsKK,
        tcnNbrOfHets=nbrOfHetsKK
      );
      verbose && cat(verbose, "Total CN segmentation table (expanded):");
      verbose && print(verbose, tcnSegmentsKK);

      # Sanity check
      stopifnot(nrow(tcnSegmentsKK) == nrow(dhSegments));

      # Combine TCN and DH segmentation results
      tcndhSegments <- cbind(
        tcnId=rep(kk, times=nrow(dhSegments)),
        dhId=seq(length=nrow(dhSegments)),
        tcnSegmentsKK,
        dhSegments
      );

      segs[[kk]] <- tcndhSegments;

      verbose && cat(verbose, "(TCN,DH) segmentation for one total CN segment:");
      verbose && print(verbose, segs[[kk]]);

      verbose && exit(verbose);
    } # for (kk ...)

    segs <- Reduce(rbind, segs);
    rownames(segs) <- NULL;
  } # if (flavor == "tcn")

  # Sanity check
  stopifnot(nrow(dhSegRows) == nrow(tcnSegsExpanded));
  rownames(tcnSegRows) <- rownames(dhSegRows) <- NULL;

  stopifnot(all(tcnSegRows[,1] <= tcnSegRows[,2], na.rm=TRUE));
  stopifnot(all(tcnSegRows[-nrow(tcnSegRows),2] < tcnSegRows[-1,1], na.rm=TRUE));
  if (flavor != "tcn") {
    stopifnot(all(dhSegRows[,1] <= dhSegRows[,2], na.rm=TRUE));
  }
  stopifnot(all(dhSegRows[-nrow(dhSegRows),2] < dhSegRows[-1,1], na.rm=TRUE));
  stopifnot(all(tcnSegsExpanded[,1] <= tcnSegsExpanded[,2], na.rm=TRUE));
##  stopifnot(all(tcnSegsExpanded[-nrow(tcnSegsExpanded),2] < tcnSegsExpanded[-1,1], na.rm=TRUE));

  # Move 'chromosome' column to the first column
  idx <- match("chromosome", names(segs));
  idxs <- c(idx, seq(length=ncol(segs))[-idx]);
  segs <- segs[,idxs,drop=FALSE];
  verbose && print(verbose, segs);

  verbose && enter(verbose, "Calculating (C1,C2) per segment");
  # Append (C1,C2) estimates
  tcn <- segs$tcnMean;
  dh <- segs$dhMean;
  C1 <- 1/2*(1-dh)*tcn;
  C2 <- tcn - C1;
  segs <- cbind(segs, c1Mean=C1, c2Mean=C2);
  verbose && exit(verbose);

  nbrOfSegs <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegs);

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create result object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- list(
    alphaTCN = alphaTCN,
    alphaDH = alphaDH,
    flavor = flavor,
    tbn = tbn,
    joinSegments = joinSegments,
    knownSegments = knownSegments,
    seed = seed
  );

  # Should we drop attributes? /HB 2010-09-24
  stopifnot(all(data$index == seq(length=nrow(data))));
  data$index <- NULL; # Drop, because it is guaranteed to be ordered
  class(data) <- c("PairedPSCNData", class(data));

  class(segs) <- c("PairedPSCNSegments", class(segs));

  fit <- list(
    data = data,
    output = segs,
    tcnSegRows = tcnSegsExpanded,
    dhSegRows = dhSegRows,
    params = params
  );

  class(fit) <- c("PairedPSCBS", "PSCBS", "AbstractCBS");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (avgTCN != "mean" || avgDH != "mean") {
    verbose && enter(verbose, "Updating mean level using different estimator");
    verbose && cat(verbose, "TCN estimator: ", avgTCN);
    verbose && cat(verbose, "DH estimator: ", avgDH);
    fit <- updateMeans(fit, avgTCN=avgTCN, avgDH=avgDH, verbose=less(verbose, 20));
    verbose && exit(verbose);
  }

  if (is.element(flavor, c("tcn&dh", "sqrt(tcn)&dh"))) {
    fit$params$flavor <- gsub("&", ",", flavor, fixed=TRUE); # AD HOC.
    fit <- postsegmentTCN(fit, verbose=verbose);

    # Sanity check
    CT <- fit$data$CT;
    tcnSegRows <- fit$tcnSegRows;
    dhSegRows <- fit$dhSegRows;
    for (jj in 1:nrow(tcnSegRows)) {
      tcnSegRowJJ <- unlist(tcnSegRows[jj,,drop=TRUE], use.names=FALSE);
      dhSegRowJJ <- unlist(dhSegRows[jj,,drop=TRUE], use.names=FALSE);
      stopifnot(
        is.na(tcnSegRowJJ[1]) || is.na(dhSegRowJJ[1]) ||
        # A TCN segment must start at or before a DH segment...
        (tcnSegRowJJ[1] <= dhSegRowJJ[1]) ||
        # ...unless there was an outlier at the left edge.
        (is.na(CT[dhSegRowJJ[1]]) && (tcnSegRowJJ[1] - 1L <= dhSegRowJJ[1]))
      );
      stopifnot(
        is.na(tcnSegRowJJ[2]) || is.na(dhSegRowJJ[2]) ||
        # A TCN segment must end at or after a DH segment...
        (dhSegRowJJ[2] <= tcnSegRowJJ[2]) ||
        # ...unless there was an outlier at the right edge.
        (is.na(CT[dhSegRowJJ[2]]) && (dhSegRowJJ[2] <= tcnSegRowJJ[2] + 1L))
      );
    } # for (jj ...)
    # Not needed anymore
    CT <- tcnSegRows <- dhSegRows <- NULL;
  }

  verbose && print(verbose, head(as.data.frame(fit)));
  verbose && print(verbose, tail(as.data.frame(fit)));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit;
}) # segmentByPairedPSCBS()



setMethodS3("segmentByPairedPSCBS", "data.frame", function(CT, ...) {
  # To please R CMD check
  data <- CT;

  segmentByPairedPSCBS(CT=data$CT, thetaT=data$thetaT, thetaN=data$thetaN,
                       betaT=data$betaT, betaN=data$betaN,
                       muN=data$muN, rho=data$rho,
                       chromosome=data$chromosome, x=data$x, ...);
})



setMethodS3("segmentByPairedPSCBS", "PairedPSCBS", function(...) {
  resegment(...);
}) # segmentByPairedPSCBS()



############################################################################
# HISTORY:
# 2014-06-08
# o Now segmentByPairedPSCBS() gives a warning about future change of the
#   default value of argument 'preserveScale' (from current TRUE to FALSE).
#   The warning only appears if the argument is not specified explicitly.
# 2014-03-30
# o As an alternative to argument 'CT', segmentByPairedPSCBS() now accepts
#   arguments 'thetaT' and 'thetaN', in case 'CT' is calculated as
#   CT=2*thetaT/thetaN (and all of 'CT', 'thetaT' and 'thetaN' are stored
#   as part the locus-level data signals.
# 2014-01-29
# o ROBUSTNESS: Now segmentByPairedPSCBS() asserts that argument 'muN'
#   is not all NAs.  Similarily, if 'muN' is called from 'betaN' the
#   same assertion is done after calling.
# 2013-02-01
# o BUG FIX: segmentByPairedPSCBS(..., avgDH="median") only worked for
#   single-chromosome data.  Same for avgTCN="median".
# 2013-01-16
# o Added arguments 'avgTCN' and 'avgDH' to segmentByPairedPSCBS().
# 2012-09-15
# o Added argument 'dropMissingCT' to segmentByPairedPSCBS().
# 2012-09-13
# o CONSISTENCY FIX: Changed the behavior of extreme values of argument
#   'undoTCN' and 'undoDH' to segmentByPairedPSCBS() such that it is
#   consistent with the new rules for 'undo' of segmentByCBS().
# 2012-07-22
# o GENERALIZATION/BUG FIX: Now segmentByPairedPSCBS() drops loci for
#   which CT is missing (regardless of betaT). For instance, in rare cases
#   when the reference (e.g. the normal) is missing, then it may be that
#   CT is missing while betaT is not.
# o Now the verbose output in segmentByPairedPSCBS() specifies the range
#   of segments with greater precision.
# 2012-04-20
# o Now it is possible to skip the DH segmentation in Paired PSCBS, i.e.
#   segmentByPairedPSCBS(..., flavor="tcn").
# o BUG FIX: segmentByPairedPSCBS() would throw "error in `$<-.data.frame
#   `(`*tmp*`, "rho" ..." if some loci points has unknown genomic positions.
# 2011-11-19
# o GENERALIZATION: Now it is possible to run Paired PSCBS (without
#   TumorBoost) when only genotypes but not BAFs are available for the
#   matched normal.
# 2011-11-17
# o ROBUSTNESS: Now all internal iterative calls to segmentByPairedPSCBS()
#   and segmentByCBS() have an explicit seed=NULL.
# o BUG FIX: Now 'tbn' argument is correctly preserved in the stored
#   parameter settings of segmentByPairedPSCBS().
# o BUG FIX: segmentByPairedPSCBS() would give an error when trying to
#   segment DH if the TCN segment contains no data points, which could
#   happen if 'knownSegments' specifies an empty segment, centromere.
# o Added segmentByPairedPSCBS() for PairedPSCBS, which is just a
#   wrapper for resegment().
# 2011-10-21
# o Now segmentByPairedCBS() handles forced separators in 'knownSegments'.
# 2011-10-20
# o CLEANUP: Dropped a stray debug output message in segmentByPairedPSCBS().
# o Replaced argument 'knownCPs' with 'knownSegments' for  segmentByCBS().
# 2011-10-02
# o Added segmentByPairedPSCBS() for data.frame such that the locus-level
#   data arguments can also be passed via a data.frame.
# 2011-09-04
# o ROBUSTNESS: Added drop=FALSE to matrix subsettings.
# o CLEANUP: Removed all references to/usage of DNAcopy fields, which
#   are no longer part of the CBS class.
# 2011-09-03
# o Updated code to not use deprecated argument 'columnNamesFlavor'
#   of segmentByCBS().
# 2011-08-08
# o BUG FIX: If dropSegmentationOutliers() would drop an outlier next to
#   a change point, such that total copy-number signal would become NA,
#   then the sanity checks that TCN segments always overlaps DH segments
#   would fail.  Now the sanity checks are aware of this special case.
# o Moved the sanity checks that tests the TCN and DH "segRows" from the
#   bootstrapTCNandDHByRegion() to segmentByPairedPSCBS().  This is the
#   first step to fix a failure in the sanity checks that could happend
#   iff one first run dropSegmentationOutliers().
# 2011-07-15
# o DOCUMENTATION: Added a section to help("segmentByPairedPSCBS") on
#   the importance of doing a whole-genome PSCBS segmentations if
#   calling AB and LOH states afterward.
# 2011-07-14
# o DOCUMENTATION: Added to the help that arguments betaT, betaN and muN
#   may contain NAs for non-polymorphic loci.
# o BUG FIX/ROBUSTNESS: In some cases, the segmentation table would
#   contain column names with incorrect capitalization, e.g. "tcnnbrOfLoci"
#   instead of "tcnNbrOfLoci".  This would cause several downstream
#   methods to give an error.  The reason for this is that the Hmisc
#   package, if loaded after R.utils, overrides capitalize() in R.utils
#   with another (buggy?) capitalize() function.  To avoid this, we
#   now everywhere specify explicitly that we want to the one in R.utils.
# 2011-07-06
# o DOCUMENTATION: The description of argument 'chromosome' for
#   segmentByPairedPSCBS() did not describe how to segment multiple
#   chromosomes in one call.
# 2011-07-05
# o BUG FIX: Output fields 'tcnNbrOfSNPs'and 'tcnNbrOfHets' were mistakenly
#   labelled as 'tcnNbrOr...'.  Thanks Christine Ho at UC Berkeley for
#   reporting on this.
# 2011-06-28
# o DOCUMENTATION: Clarified that argument 'CT' should be tumor copy
#   number ratios relative to the normal.
# 2011-06-14
# o CONVENTION: Changed the column names of returned data frames.
#   They now follow the camelCase naming convention and are shorter.
# 2011-05-29
# o Renamed options to reflect new package name.
# 2011-04-07
# o ROBUSTNESS: Added validation of the different 'tcnSegRows' and
#   'dhSegRows' calculations in segmentByPairedPSCBS().  This helps
#   us narrow down a bug in postsegmentTCN().
# 2010-12-09
# o BUG FIX: When there were multiple chromsomes processed by
#   segmentByPairedPSCBS(), then the returned data object would
#   contain 'betaT' identical to 'betaTN'.
# 2010-12-02
# o Now segmentByPairedPSCBS() uses option "psCBS/sanityChecks/tolerance".
# 2010-11-30
# o Now segmentByPairedPSCBS() returns data frames 'tcnLociToExclude'
#   and 'dhLociToExclude'.
# o BUG FIX: Argument 'flavor' of segmentByPairedPSCBS() would be ignored
#   if multiple chromsomomes were segmented.
# 2010-11-28
# o BUG FIX: Iff argument 'chromosome' to segmentByPairedPSCBS() was of
#   length greater than one and specified exactly one unique chromosome,
#   then exception "Number of elements in argument 'chromosome' should
#   be exactly 8712 not 86209 value(s)" would be thrown.
# 2010-11-27
# o BUG FIX: segmentByPairedPSCBS() would not accept missing values in
#   argument 'chromosome'.
# o Now arguments '...' of segmentByPairedPSCBS() are passed to
#   the two segmentByCBS() calls.
# 2010-11-22
# o BUG FIX: segmentByPairedPSCBS() would not subset the correct set of
#   DH signals if there were some missing values in TCN.
# 2010-11-21
# o Changed the default to flavor="tch&dh".
# o Added support for flavors "tcn&dh", which, contrary to "tcn,dh",
#   enforces TCN and DH to have the same change points.
# o Now segmentByPairedPSCBS() also returns minor and major copy numbers
#   for each segment.
# o Forgot to return arguments 'joinSegments' & 'knownCPs' in 'params'.
# 2010-11-20
# o Now it is possible to specify the boundaries of the regions to be
#   segmented as known change points via argument 'knownCPs'.
# o Added argument 'joinSegments' to segmentByPairedPSCBS() in order to
#   specify if neighboring segments should be joined or not.
# o Now segmentByCBS() allows for unknown genomic positions.
# o Now segmentByCBS() allows also for missing total CN signals.
# 2010-11-16
# o BUG FIX: In the rare cases where two loci at the same positions are
#   split up into two neighboring segments, then segmentByPairedPSCBS()
#   would fail to infer which they were if and only if the loci were not
#   ordered along the genome.  This could happen with for instance
#   Affymetrix GenomeWideSNP_6 data.
# o DOCUMENTATION: Clarified the form of argument 'muN', and added
#   references to papers and cross links to more internal methods.
# 2010-11-04
# o BUG FIX: There was a stray/debug stop() statement left in
#   segmentByPairedPSCBS() causing an "error" in the rare case
#   when loci that have the same physical locations are split
#   into two different segments.
# 2010-11-02
# o Added arguments 'undoTCN' and 'undoDH' to segmentByPairedPSCBS().
# o BUG FIX: Arguments 'alphaTCN' and 'alphaDH' of segmentByPairedPSCBS()
#   were not used when more than one chromosome were segmented.
# 2010-10-25
# o BUG FIX: Now the correct set of loci are extracted from each TCN
#   segment, in the rare case that two neighboring TCN segments have
#   the same end points.
# 2010-10-18
# o Added arguments 'alphaTCN' and 'alphaDH' to segmentByPairedPSCBS().
# o Now segmentByPairedPSCBS() can segment multiple chromosomes.
# 2010-10-17
# o Added argument 'tbn' to segmentByPairedPSCBS() specifying whether
#   TumorBoostNormalization should be applied or not.
# 2010-10-10
# o The default is now to segment TCN on the original scale, not the sqrt().
# o Added flavor "sqrt(tcn),dh", which is segments sqrt(TCN) and then DH,
#   as original proposed by ABO.
# 2010-10-03
# o CLEAN UP: Now segmentByPairedPSCBS() is making use of argument
#   'chromosome' of segmentByCBS().
# 2010-10-02
# o Argument 'chromosome' default to 0 and have to be a finite integer.
# 2010-09-24
# o Now the 'data' field returned is a data.frame (no longer a list).
# o Now the 'chromosome' field of the data field is expanded to have the
#   same number of elements as the other locus fields.
# 2010-09-18
# o Added argument 'chromosome' to segmentByPairedPSCBS(), which, if given,
#   adds a chromosome column to the data and segmentation results.
# 2010-09-08
# o Now segmentByPairedPSCBS() also returns the TumorBoost normalized data.
#   This also means that plot() for PairedPSCBS no longer has to
#   recalculate them.
# 2010-09-04
# o Added drawLevels() for PairedPSCBS.
# o Added as.data.frame() and print() for PairedPSCBS.
# 2010-09-03
# o Added plot().
# 2010-07-09
# o The segmentByPairedPSCBS() method was written completely from scratch.
# o Created.
############################################################################
