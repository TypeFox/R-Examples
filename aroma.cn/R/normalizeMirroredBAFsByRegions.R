###########################################################################/**
# @set "class=matrix"
# @RdocMethod normalizeMirroredBAFsByRegions
# @alias normalizeMirroredBAFsByRegions
#
# @title "Normalizes region-level mirrored allele B fractions (mBAFs)"
#
# \description{
#  @get "title" for heterozygous and homozygous SNPs by rescaling both
#  equally much such that the normalized homozygous mBAFs become exactly
#  or close to their expected values (which is at zero and one).
# }
#
# @synopsis
#
# \arguments{
#   \item{data}{A @numeric Kx2 or Kx3 @matrix, where K is the number of
#     segments and the first and the second column contains average
#     heterozygous and homozygous mBAF estimates, respectively.
#     The third column, which is optional, contains total copy numbers.}
#   \item{flavor}{A @character string specifying how the normalization
#     function/scale factors are estimated.}
#   \item{...}{Additional arguments passed @see "aroma.light::fitXYCurve",
#     which is used if \code{flavor="total"}.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   A @numeric @matrix of the same dimensions as argument \code{data}.
# }
#
# \details{
#   The mBAFs for heterozygous SNPs are also known as the
#   Decrease in Heterozygosity signals (DHs).
# }
#
# %examples "../incl/normalizeMirroredBAFsByRegions.Rex"
#
# @author "HB, PN"
#
# @keyword internal
#*/###########################################################################
setMethodS3("normalizeMirroredBAFsByRegions", "matrix", function(data, flavor=c("plain", "total"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'data':
  dim <- dim(data);
  if (!is.element(dim[2], c(2,3))) {
    dimStr <- sprintf("%dx%d", dim[1], dim[2]);
    throw("Argument 'data' must be an Kx2 or Kx3 matrix: ", dimStr);
  }
  nbrOfSegments <- nrow(data);

  # Argument 'flavor':
  flavor <- match.arg(flavor);
  if (flavor == "total") {
    if (dim[2] < 3) {
      throw("Argument 'data' must contain 3 column when flavor=\"total\": ", dim[2]);
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Rescaling segment-level mirrored allele B fractions (mBAFs)");

  verbose && cat(verbose, "Flavor: ", flavor);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);
  verbose && cat(verbose, "Mirrored allele B fractions (data):");
  verbose && str(verbose, data);
  verbose && summary(verbose, data);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate normalization function/scale factors
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  modelFit <- list(scale=NA);
  if (flavor == "plain") {
    verbose && enter(verbose, "Estimating scale factors for each segment independently");

    # Estimate the scale factors for each segment independently
    y <- data[,2,drop=TRUE];
    scale <- 1 / y;
    verbose && cat(verbose, "Scale factors (one per segment):");
    verbose && str(verbose, scale);
    # Not needed anymore
    y <- NULL;

    verbose && exit(verbose);
  } else if (flavor == "total") {
    verbose && enter(verbose, "Estimating scale factors globally as a function of total copy numbers");

    verbose && enter(verbose, "Fit global relationship between homozygous mBAFs and TCNs");
    X <- data[,c(3,2),drop=FALSE];
    verbose && cat(verbose, "(TCN,mBAFhom):");
    verbose && str(verbose, X);
    fit <- .fitXYCurve(X, ..., verbose=verbose);
    # Not needed anymore
    X <- NULL;
    verbose && str(verbose, fit);
    verbose && exit(verbose);

    verbose && enter(verbose, "Predict scale factors");
    yHat <- fit$predictY(data[,3,drop=TRUE]);
    y <- data[,2];
    plot(y,yHat, xlim=c(0,1), ylim=c(0,1));
    scale <- 1 / yHat;
    verbose && cat(verbose, "Predicted scale factors (one per segment):");
    verbose && str(verbose, scale);
    verbose && exit(verbose);

    # Store the normalization function
    modelFit$subFit <- fit;

    # Not needed anymore
    fit <- yHat <- NULL;

    verbose && exit(verbose);
  }

  # Store the scale estimates
  modelFit$scale <- scale;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalizing mBAFs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing the mirrored allele B fractions");
  dataN <- data;
  dataN[,1:2] <- scale * data[,1:2];
  verbose && str(verbose, dataN);
  verbose && summary(verbose, dataN);
  verbose && exit(verbose);

  attr(dataN, "modelFit") <- modelFit;

  verbose && exit(verbose);

  # Sanity check
  stopifnot(dim(dataN) == dim);

  dataN;
}) # normalizeMirroredBAFsByRegions()


##############################################################################
# HISTORY
# 2012-04-16
# o normalizeMirroredBAFsByRegions() now explicitly require 'aroma.light'.
# 2010-09-08 [PN+HB]
# o Added normalizeMirroredBAFsByRegions().
# o Created.
##############################################################################
