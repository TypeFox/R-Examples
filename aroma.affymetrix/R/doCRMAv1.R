###########################################################################/**
# @RdocDefault doCRMAv1
# @alias doCRMAv1.AffymetrixCelSet
# @alias doASCRMAv1
# @alias doASCRMAv1.default
#
# @title "Estimation and assessment of raw copy numbers at the single locus level (CRMA v1)"
#
# \description{
#  @get "title" based on [1].
#  The algorithm is processed in bounded memory, meaning virtually
#  any number of arrays can be analyzed on also very limited computer
#  systems.
# }
#
# \usage{
#   @usage doCRMAv1,AffymetrixCelSet
#   @usage doCRMAv1,default
#   @usage doASCRMAv1,default
# }
#
# \arguments{
#  \item{csR, dataSet}{An @see "AffymetrixCelSet" (or the name of an @see "AffymetrixCelSet").}
#  \item{shift}{An tuning parameter specifying how much to shift the
#   probe signals before probe summarization.}
#  \item{combineAlleles}{A @logical specifying whether allele probe pairs
#   should be summed before modelling or not.}
#  \item{lengthRange}{An optional @numeric vector of length two passed
#   to @see "FragmentLengthNormalization".}
#  \item{arrays}{A @integer @vector specifying the subset of arrays
#   to process.  If @NULL, all arrays are considered.}
#  \item{drop}{If @TRUE, the summaries are returned, otherwise
#   a named @list of all intermediate and final results.}
#  \item{verbose}{See @see "Verbose".}
#  \item{...}{Additional arguments used to set up @see "AffymetrixCelSet" (when argument \code{dataSet} is specified).}
# }
#
# \value{
#   Returns a named @list, iff \code{drop == FALSE}, otherwise
#   only @see "ChipEffectSet" object.
# }
#
#
# \section{Allele-specific or only total-SNP signals}{
#   If you wish to obtain allele-specific estimates for SNPs, which
#   are needed to call genotypes or infer parent-specific copy numbers,
#   then use argument \code{combineAlleles=FALSE}.  Total copy number
#   signals are still available.
#   If you know for certain that you will not use allele-specific
#   estimates, you will get slightly less noisy signals
#   (very small difference) if you use \code{combineAlleles=TRUE}.
#
#   \code{doASCRMAv1(...)} is a wrapper for
#   \code{doCRMAv1(..., combineAlleles=FALSE)}.
# }
#
# \references{
#  [1] H. Bengtsson, R. Irizarry, B. Carvalho & T.P. Speed.
#      \emph{Estimation and assessment of raw copy numbers at the
#      single locus level},
#      Bioinformatics, 2008.\cr
# }
#
# \seealso{
#  For CRMA v2 (recommended by authors), which is a single-array
#  improvement over CRMA v1, see @see "doCRMAv2".
# }
#
# @author "HB"
#*/###########################################################################
setMethodS3("doCRMAv1", "AffymetrixCelSet", function(csR, shift=+300, combineAlleles=TRUE, lengthRange=NULL, arrays=NULL, drop=TRUE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csR':
  className <- "AffymetrixCelSet";
  if (!inherits(csR, className)) {
    throw(sprintf("Argument 'csR' is not a %s: %s", className, class(csR)[1]));
  }

  # Argument 'shift':
  shift <- Arguments$getNumeric(shift);

  # Argument 'combineAlleles':
  combineAlleles <- Arguments$getLogical(combineAlleles);

  # Argument 'arrays':
  if (!is.null(arrays)) {
    throw("Not supported. Argument 'arrays' should be NULL.");
    arrays <- Arguments$getIndices(arrays, max=length(csR));
  }

  # Argument 'drop':
  drop <- Arguments$getLogical(drop);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "CRMAv1");
  verbose && cat(verbose, "Arguments:");
  verbose && cat(verbose, "combineAlleles: ", combineAlleles);
  arraysTag <- seqToHumanReadable(arrays);
  verbose && cat(verbose, "arrays:");
  verbose && str(verbose, arraysTag);

  # Backward compatibility
  ram <- list(...)$ram;
  if (!is.null(ram)) {
    ram <- Arguments$getDouble(ram, range=c(0,Inf));
    verbose && cat(verbose, "ram: ", ram);
    warning("Argument 'ram' of doCRMAv1() is deprecated. Instead use setOption(aromaSettings, \"memory/ram\", ram).");
    oram <- setOption(aromaSettings, "memory/ram", ram);
    on.exit({
      setOption(aromaSettings, "memory/ram", oram);
    });
  }

  # List of objects to be returned
  res <- list();
  if (!drop) {
    res <- c(res, list(csR=csR));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup data set to be processed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Data set");
  verbose && print(verbose, csR);

  if (!is.null(arrays)) {
    verbose && enter(verbose, "CRMAv1/Extracting subset of arrays");
    csR <- extract(csR, arrays, onDuplicates="error");
    verbose && cat(verbose, "Data subset");
    verbose && print(verbose, csR);
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # CRMAv1
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "CRMAv1/Allelic crosstalk calibration");
  acc <- AllelicCrosstalkCalibration(csR, model="CRMA", tags="*,v1");
  verbose && print(verbose, acc);
  csC <- process(acc, verbose=verbose);
  verbose && print(verbose, csC);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(acc=acc, csC=csC));
  }

  # Clean up
  # Not needed anymore
  csR <- acc <- NULL;
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "CRMAv1/Probe summarization");
  plm <- RmaCnPlm(csC, mergeStrands=TRUE, combineAlleles=combineAlleles,
                                                            shift=shift);
  verbose && print(verbose, plm);
  if (length(findUnitsTodo(plm)) > 0) {
    # Fit CN probes quickly (~5-10s/array + some overhead)
    units <- fitCnProbes(plm, verbose=verbose);
    verbose && str(verbose, units);
    # Fit remaining units, i.e. SNPs (~5-10min/array)
    units <- fit(plm, verbose=verbose);
    verbose && str(verbose, units);
    # Not needed anymore
    units <- NULL;
  }
  verbose && print(verbose, gc);
  ces <- getChipEffectSet(plm);
  verbose && print(verbose, ces);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(ces=ces, plm=plm));
  }

  # Clean up
  # Not needed anymore
  plm <- csC <- NULL;
  gc <- gc();

  verbose && enter(verbose, "CRMAv1/PCR fragment-length normalization");
  fln <- FragmentLengthNormalization(ces, target="zero", lengthRange=lengthRange);
  verbose && print(verbose, fln);
  cesN <- process(fln, verbose=verbose);
  verbose && print(verbose, cesN);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(fln=fln, cesN=cesN));
  }

  # Clean up
  # Not needed anymore
  fln <- ces <- NULL;
  gc <- gc();

  verbose && enter(verbose, "CRMAv1/Export to technology-independent data files");
  dsNList <- exportTotalAndFracB(cesN, verbose=verbose);
  verbose && print(verbose, dsNList);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(dsNList=dsNList));
  }

  # Clean up
  # Not needed anymore
  cesN <- NULL;
  gc <- gc();

  verbose && exit(verbose);

  # Return only the final results?
  if (drop) {
    res <- dsNList;
  }

  res;
}) # doCRMAv1()


setMethodS3("doCRMAv1", "default", function(dataSet, ..., verbose=FALSE) {
  .require <- require
  .require("aroma.affymetrix") || throw("Package not loaded: aroma.affymetrix")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "CRMAv1");

  verbose && enter(verbose, "CRMAv1/Setting up CEL set");
  csR <- AffymetrixCelSet$byName(dataSet, ..., verbose=less(verbose, 50),
                                                  .onUnknownArgs="ignore");
  verbose && print(verbose, csR);
  verbose && exit(verbose);

  dsNList <- doCRMAv1(csR, ..., verbose=verbose);

  # Clean up
  # Not needed anymore
  csR <- NULL;
  gc <- gc();

  verbose && exit(verbose);

  dsNList;
})


setMethodS3("doASCRMAv1", "default", function(...) {
  .require <- require
  .require("aroma.affymetrix") || throw("Package not loaded: aroma.affymetrix")

  doCRMAv1(..., combineAlleles=FALSE)
})


############################################################################
# HISTORY:
# 2013-05-02
# o Removed argument 'ram' in favor of aroma option 'memory/ram'.
# 2012-09-05
# o ROBUSTNESS: Now doCRMAv1() adds also tag "v1" to the allele-specific
#   calibration step.  The reason for this is to differentiate it from
#   the output of doCRMAv2().  NOTE: This update means that any old CRMAv1
#   analyses will not be detected by doCRMAv1(); to have doCRMAv1() detect
#   those add tag "v1" in that calibration step, e.g. "ACC,-XY,v1".
# 2011-04-07
# o Added argument 'drop'.
# 2011-03-14
# o doCRMAv1() gained argument 'lengthRange', which is passed to
#   the constructor of FragmentLengthNormalization.
# 2010-06-21
# o Added doASCRMAv1() for a convenient allele-specific CRMAv1 wrapper.
# 2010-06-07
# o Added argument shift=+300 to doCRMAv1().
# 2010-05-17
# o CORRECTION: doCRMAv1() forgot to shift +300 the signals before
#   doing the probe-level summarization.
# 2010-04-21
# o BUG FIX: doCRMAv1() for AffymetrixCelSet used undefined 'csN' internally
#   instead of 'csC'.
# 2010-04-04
# o Created from doCRMAv2.R.
# o (Re)created.
############################################################################
