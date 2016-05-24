###########################################################################/**
# @set "class=array"
# @RdocMethod calmateByThetaAB
# @alias calmateByThetaAB
# 
# @title "Normalize allele-specific copy numbers (CA,CB)"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{data}{An Jx2xI @numeric @array, where J is the number of SNPs,
#          2 is the number of alleles, and I is the number of samples.}
#  \item{references}{An index @vector in [1,I] or a @logical @vector 
#     of length I specifying which samples are used when calculating the
#     reference signals.  If @NULL, all samples are used. At least 3 samples.}
#  \item{...}{Additional arguments passed to the internal fit function
#     @see "fitCalMaTeInternal".}
#  \item{truncate}{If @TRUE, final ASCNs are forced to be non-negative
#     while preserving the total CNs.}
#  \item{refAvgFcn}{(optional) A @function that takes a JxI @numeric @matrix
#     an argument \code{na.rm} and returns a @numeric @vector of length J.
#     It should calculate some type of average for each of the J rows, e.g.
#     @see "matrixStats::rowMedians".  
#     If specified, then the total copy numbers of the calibrated ASCNs
#     are standardized toward (twice) the average of the total copy numbers
#     of the calibrated reference ASCNs.}
#  \item{flavor}{A @character string specifying which flavor of the
#     CalMaTe algorithm to use for fitting the model.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an Jx2xI @numeric @array
#   with the same dimension names as argument \code{data}.
# }
#
# \section{Flavors}{
#   For backward compatibility, we try to keep all major versions of
#   the CalMaTe algorithm available.  Older versions can be used by
#   specifying argument \code{flavor}.
#   The default flavor is \code{v2}.
#   For more information about the different flavors, 
#   see @see "fitCalMaTeInternal".
# }
#
# @examples "../incl/calmateByThetaAB.Rex"
#
# \references{
#  [1] @include "../incl/OrtizM_etal_2012.Rd" \cr 
# }
#
# \seealso{
#  To calibrate (total,fracB) data, 
#  see @seemethod "calmateByTotalAndFracB".
#  We strongly recommend to always work with (total,fracB) data
#  instead of (CA,CB) data, because it is much more general.
#
#  For further information on the internal fit functions, see
#  @see "fitCalMaTeInternal".
# }
#*/###########################################################################
setMethodS3("calmateByThetaAB", "array", function(data, references=NULL, ..., truncate=FALSE, refAvgFcn=NULL, flavor=c("v2", "v1"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'data':
  if (!is.array(data)) {
    throw("Argument 'data' is not an array: ", class(data)[1]);
  }
  dim <- dim(data);
  dimnames <- dimnames(data);
  if (length(dim) != 3) {
    throw("Argument 'data' is not a 3-dimensional array: ", 
                                                paste(dim, collapse="x"));
  }
  if (dim[2] != 2) {
    throw("Argument 'data' is not a Jx2xI-dimensional array: ", 
                                                paste(dim, collapse="x"));
  }
  if (!is.null(dimnames[[2]])) {
    if (!identical(dimnames[[2]], c("A", "B"))) {
      throw("If given, the names of the allele (2nd) dimension of the Jx2xI-dimensional array (argument 'data') have to be 'A' & 'B': ", paste(dimnames[[2]], collapse=", "));
    }
  }

  nbrOfSamples <- dim[3];
  if (nbrOfSamples < 3) {
    throw("Argument 'data' contains less than three samples: ", nbrOfSamples);
  }

  # Argument 'references':
  if (is.null(references)) {
    # The default is that all samples are used to calculate the reference.
    references <- seq(length=nbrOfSamples);
  } else if (is.logical(references)) {
    if (length(references) != nbrOfSamples) {
      throw("Length of argument 'references' does not match the number of samples in argument 'data': ", length(references), " != ", nbrOfSamples);
    }
    references <- which(references);
  } else if (is.numeric(references)) {
    references <- as.integer(references);
    if (any(references < 1 | references > nbrOfSamples)) {
      throw(sprintf("Argument 'references' is out of range [1,%d]: %d", nbrOfSamples), length(references));
    }
  }
  if (length(references) < 3) {
    throw("Argument 'reference' specify less than three reference samples: ", length(references));
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);    



  # From here on we force dimension names on the 2nd dimension
  dimnames(data)[[2]] <- c("A", "B");

  
  verbose && enter(verbose, "calmateByThetaAB()");
  verbose && cat(verbose, "ASCN signals:");
  verbose && str(verbose, data);
  verbose && cat(verbose, "Reference samples:");
  verbose && str(verbose, references);

  verbose && enter(verbose, "Identifying non-finite data points");
  # Keep finite values
  ok <- (is.finite(data[,"A",,drop=FALSE]) & is.finite(data[,"B",,drop=FALSE]));
  dim(ok) <- dim(ok)[-2]; # Drop 2nd dimension
  ok <- rowAlls(ok);
  verbose && summary(verbose, ok);
  hasNonFinite <- any(!ok);
  if (hasNonFinite) {
    verbose && enter(verbose, "Excluding non-finite data points");
    dataS <- data[ok,,,drop=FALSE];
    verbose && str(verbose, data);
    verbose && exit(verbose);
    dim <- dim(dataS);
  } else {
    verbose && cat(verbose, "All data points are finite.");
    dataS <- data;
  }
  verbose && exit(verbose);


  verbose && enter(verbose, "Fitting CalMaTe");
  verbose && cat(verbose, "Algorithm flavor: ", flavor);
  if (flavor == "v2") {
    fitFcn <- fitCalMaTeV2;
  } else if (flavor == "v1") {
    fitFcn <- fitCalMaTeV1;
  } else {
    throw("Unknown algorithm flavor: ", flavor);
  }

  nbrOfSNPs <- dim(dataS)[1];
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSNPs);
  verbose && printf(verbose, "Number of SNPs left: ");
  # Drop dimnames for faster processing
  dimnames(dataS) <- NULL;
  # Used for sanity check inside loop
  dimCjj <- dim(dataS)[-1];
  for (jj in seq(length=nbrOfSNPs)) {
    if (verbose && (jj %% 500 == 1)) {
      writeRaw(verbose, sprintf("%d, ", nbrOfSNPs-jj+1));
    }
    Cjj <- dataS[jj,,,drop=FALSE];  # An 1x2xI array
    dim(Cjj) <- dimCjj;             # A 2xI matrix
    CCjj <- fitFcn(Cjj, references=references, ...);
    # Sanity check
    stopifnot(identical(dim(CCjj), dimCjj));
    dataS[jj,,] <- CCjj;
  } # for (jj ...)
  if (verbose) writeRaw(verbose, "done.\n");
  verbose && exit(verbose);

  if (hasNonFinite) {
    verbose && enter(verbose, "Expanding to array with non-finite");
    dataC <- data;
    dataC[ok,,] <- dataS;
    verbose && str(verbose, dataC);
    verbose && exit(verbose);
  } else {
    dataC <- dataS;
    dimnames(dataC) <- dimnames(data);
  }
  rm(dataS);

  # Sanity check
  stopifnot(identical(dim(dataC), dim(data)));

  verbose && cat(verbose, "Calibrated ASCN signals:");
  verbose && str(verbose, dataC);

  if (truncate){
    dataC <- truncateThetaAB(dataC);
    verbose && cat(verbose, "Truncated ASCN signals:");
    verbose && str(verbose, dataC);
  } 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Standardize toward a custom average of the references?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(refAvgFcn)) {
    verbose && enter(verbose, "Standardize total copy numbers toward the average reference signals");
    # Extract reference signals
    dataCR <- dataC[,,references,drop=FALSE];
    # Calculate total copy number signals
    yCR <- dataCR[,1,,drop=FALSE]+dataCR[,2,,drop=FALSE];
    dim(yCR) <- dim(yCR)[-2]; # Drop 2nd dimension
    # Calculate the average
    yCR <- refAvgFcn(yCR, na.rm=TRUE);
    # Standardize ASCNs to this average
    dataC <- 2 * dataC / yCR;
    verbose && exit(verbose);
  }

  # Enforce the same dimension names as the input data
  dimnames(dataC) <- dimnames;

  verbose && cat(verbose, "Calibrated (A,B) signals:");
  verbose && str(verbose, dataC);

  verbose && exit(verbose);

  dataC;
}) # calmateByThetaAB()


###########################################################################
# HISTORY:
# 2012-02-20 [HB]
# o Minor speed up of calmateByThetaAB() by precalculating dim(Cjj)
#   outside the loop.
# o Made the progress verbose messages of calmateByThetaAB() tighter.
# 2012-02-19 [HB]
# o BACKWARD COMPATIBILITY: Added argument 'flavor' to calmateByThetaAB()
#   to be able to use previous versions of CalMaTe model estimators.
#   Flavor "v2" was introduced 2011-12-05.
# o SPEEDUP: the internal CalMaTe fit function is now called directly,
#   which avoids the method dispatch overhead that otherwise applies 
#   to each SNP fitted.
# 2011-12-15 [HB]
# o CLEANUP: Tidied up the validation of argument 'references' and
#   improved the corresponding error messages.
# 2011-12-07 [MO]
# o At least 3 reference samples.
# 2011-03-18 [HB]
# o BUG FIX: calmateByThetaAB() required that the 2nd dimension
#   of argument 'data' had names "A" and "B".
# 2010-08-05 [HB]
# o ROBUSTNESS: Now calmateByThetaAB() asserts that there is at least
#   two samples.
# o BUG FIX: calmateByThetaAB() would not work with only one unit or only
#   one sample.
# 2010-08-02 [HB]
# o Added argument 'refAvgFcn' to calmateByThetaAB().
# 2010-06-19 [HB]
# o Now calmateByThetaAB() uses internal truncateThetaAB() to truncate
#   (CA,CB) values.  Since it operates on arrays, it is much faster.
# 2010-06-18 [HB]
# o Now calmateByThetaAB() calls internal fitCalMaTe().
# o Added argument 'references' to calmateByThetaAB().
# 2010-06-04 [MO]
# o Created.
###########################################################################
