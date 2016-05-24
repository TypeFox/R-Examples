###########################################################################/**
# @RdocDefault segmentByNonPairedPSCBS
# @alias segmentByNonPairedPSCBS.data.frame
# @alias segmentByNonPairedPSCBS.PairedPSCBS
# @alias segmentByNonPairedPSCBS
#
# @title "Segment total copy numbers and allele B fractions using the Non-paired PSCBS method"
#
# \description{
#  @get "title" [1].
#  This method does not requires matched normals.
#  This is a low-level segmentation method.
#  It is intended to be applied to one tumor sample at the time.
# }
#
# @synopsis
#
# \arguments{
#   \item{CT}{A @numeric @vector of J tumor total copy number (TCN)
#        ratios in [0,+@Inf) (due to noise, small negative values are
#        also allowed).  The TCN ratios are typically scaled such that
#        copy-neutral diploid loci have a mean of two.}
#   \item{betaT}{A @numeric @vector of J tumor allele B fractions (BAFs)
#        in [0,1] (due to noise, values may be slightly outside as well)
#        or @NA for non-polymorphic loci.}
#   \item{...}{Additional arguments passed to @see "segmentByPairedPSCBS".}
#   \item{flavor}{A @character specifying what type of segmentation and
#     calling algorithm to be used.}
#   \item{tauA, tauB}{Lower and upper thresholds (\code{tauA < tauB} for
#     calling SNPs  heterozygous based on the tumor allele B fractions
#     (\code{betaT}).  If @NA, then they are estimates from data.
#   }
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the segmentation results as a @see "NonPairedPSCBS" object.
# }
#
# \details{
#   Internally @see "segmentByPairedPSCBS" is used for segmentation.
#   This segmentation method does \emph{not} support weights.
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
# \section{Non-Paired PSCBS with known genotypes}{
#   If allele B fractions for the matched normal (\code{betaN}) are
#   not available, but genotypes (\code{muN}) are, then it is possible
#   to run Paired PSCBS.   See @see "segmentByPairedPSCBS" for details.
# }
#
# @examples "../incl/segmentByNonPairedPSCBS.Rex"
#
# @author "HB"
#
# \references{
#  [1] @include "../incl/OlshenA_etal_2011.Rd" \cr
#  [2] @include "../incl/BengtssonH_etal_2010.Rd" \cr
# }
#
# \seealso{
#   To segment paired tumor-normal total copy numbers and allele B fractions,
#   see @see "segmentByPairedPSCBS".
#
#   To segment total copy numbers, or any other unimodal signals,
#   see @see "segmentByCBS".
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("segmentByNonPairedPSCBS", "default", function(CT, betaT, ..., flavor=c("tcn", "tcn&dh", "tcn,dh", "sqrt(tcn),dh", "sqrt(tcn)&dh"), tauA=NA, tauB=1-tauA, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'CT':
  disallow <- c("Inf");
  CT <- Arguments$getDoubles(CT, disallow=disallow);
  nbrOfLoci <- length(CT);
  length2 <- rep(nbrOfLoci, times=2);

  # Argument 'betaT':
  betaT <- Arguments$getDoubles(betaT, length=length2, disallow="Inf");

  # Argument 'flavor':
  flavor <- match.arg(flavor);
  knownFlavors <- eval(formals(segmentByPairedPSCBS.default)$flavor);
  if (!is.element(flavor, knownFlavors)) {
    throw("Segmentation flavor is not among the supported ones (", paste(sprintf("\"%s\"", knownFlavors), collapse=", "), "): ", flavor);
  }

  # Argument 'tauA' & 'tauB':
  if (!is.na(tauA) && !is.na(tauB)) {
    tauA <- Arguments$getDouble(tauA);
    tauB <- Arguments$getDouble(tauB);
    if (tauB < tauA) {
      throw("Argument 'tauA' must be smaller than 'tauB': ", tauA, " > ", tauB);
    }
    tauA <- Arguments$getDouble(tauA, range=c(-0.5, +0.5));
    tauB <- Arguments$getDouble(tauB, range=c(+0.5, +1.5));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Segmenting non-paired tumor signals using Non-paired PSCBS");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup input data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);

  # SNPs are identifies as those loci that have non-missing 'betaT'
  isSnp <- !is.na(betaT);
  nbrOfSnps <- sum(isSnp);
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSnps);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call tumor "genotypes"
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calling \"genotypes\" from tumor allele B fractions");
  verbose && str(verbose, betaT);

  if (is.na(tauA) && is.na(tauB)) {
    mBAF <- abs(betaT - 1/2);
    findPeaksAndValleys <- .use("findPeaksAndValleys", package="aroma.light");
    fitT <- findPeaksAndValleys(mBAF);
    type <- NULL; rm(list="type"); # To please 'R CMD check'.
    fitT <- subset(fitT, type == "peak");
    o <- order(fitT$density, decreasing=TRUE);
    fitT <- fitT[o,];
    fitT <- fitT[1,];
    z <- mBAF[mBAF >= fitT$x] - fitT$x;
    q <- quantile(z, probs=0.95, na.rm=TRUE, names=FALSE);
    qU <- fitT$x+q;
    verbose && cat(verbose, "Upper quantile: ", qU);
    qL <- fitT$x - q;
    verbose && cat(verbose, "Symmetric lower quantile: ", qL);
    tauA <- 1/2-qL;
    tauB <- 1/2+qL;
    verbose && cat(verbose, "(tauA, tauB) estimates: (%g,%g)", tauA, tauB);

    # Sanity check on (tauA, tauB) estimates
    if (tauB < tauA) {
      throw("Failed to estimate (tauA, tauB). The estimate 'tauA' is greater than 'tauB', which it should not: ", tauA, " > ", tauB);
    }
    tauA <- Arguments$getDouble(tauA, range=c(-0.5, +0.5));
    tauB <- Arguments$getDouble(tauB, range=c(+0.5, +1.5));
  }

  verbose && cat(verbose, "Homozygous treshholds:");
  verbose && print(verbose, c(tauA, tauB));

  isHomA <- isSnp & (betaT <= tauA);
  isHomB <- isSnp & (betaT >= tauB);
  isHom <- (isHomA | isHomB);
  isHet <- isSnp & !isHom;

  # Tumor proxy for germline genotypes
  naValue <- NA_real_;
  muNx <- rep(naValue, times=length(betaT));
  muNx[isHomA] <-   0;
  muNx[isHet]  <- 1/2;
  muNx[isHomB] <-   1;
  # Not needed anymore
  isHomA <- isHomB <- isHom <- isHet <- NULL;

  verbose && cat(verbose, "Inferred germline genotypes (via tumor):");
  verbose && str(verbose, muNx);
  verbose && print(verbose, table(muNx));
  verbose && exit(verbose);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segment using Paired PSCBS segmentation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Segment using Paired PSCBS");
  fit <- segmentByPairedPSCBS(CT=CT, betaT=betaT, muN=muNx, tbn=FALSE, flavor=flavor, ..., verbose=verbose);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Coerce fit object to Non-Paired PSCBS results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Coercing to Non-Paired PSCBS results");

  data <- fit$data;
  class(data) <- gsub("PairedPSCNData", "NonPairedPSCNData", class(data), fixed=TRUE);
#  class(data) <- c("NonPairedPSCNData", class(data));
  fit$data <- data;
  # Not needed anymore
  data <- NULL;

  segs <- fit$output;
  class(segs) <- gsub("PairedPSCNSegments", "NonPairedPSCNSegments", class(segs), fixed=TRUE);
#  class(segs) <- c("NonPairedPSCNSegments", class(segs));
  fit$output <- segs;
  # Not needed anymore
  segs <- NULL;

  params <- fit$params;
  params$tauA <- tauA;
  params$tauB <- tauB;
  fit$params <- params;
  # Not needed anymore
  params <- NULL;

#  class(fit) <- gsub("PairedPSCBS", "NonPairedPSCBS", class(fit), fixed=TRUE);
  class(fit) <- c("NonPairedPSCBS", class(fit));

  verbose && exit(verbose);


  verbose && print(verbose, head(as.data.frame(fit)));
  verbose && print(verbose, tail(as.data.frame(fit)));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit;
}) # segmentByNonPairedPSCBS()



setMethodS3("segmentByNonPairedPSCBS", "data.frame", function(CT, ...) {
  # To please R CMD check
  data <- CT;


  segmentByNonPairedPSCBS(CT=data$CT, betaT=data$betaT,
                          chromosome=data$chromosome, x=data$x, ...);
})



setMethodS3("segmentByNonPairedPSCBS", "PairedPSCBS", function(...) {
  resegment(...);
})



############################################################################
# HISTORY:
# 2013-07-19
# o ROBUSTNESS: Added a sanity check on the estimates of (tauA, tauB)
#   when they are estimated from data in segmentByNonPairedPSCBS().
# 2012-11-05
# o DOCUMENTATION FIX: example(segmentByNonPairedPSCBS) was for the
#   paired case.
# 2012-08-20
# o BUG FIX: segmentByNonPairedPSCBS() forgot to specify namespace
#   aroma.light when trying to call findPeaksAndValleys().
# 2012-04-20
# o Created from segmentByPairedPSCBS.R.
############################################################################
