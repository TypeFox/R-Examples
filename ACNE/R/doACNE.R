###########################################################################/**
# @RdocDefault doACNE
# @alias doACNE.AffymetrixCelSet
#
# @title "(ACNE)"
#
# \description{
#  @get "title" based on [1].
#  The algorithm is processed in bounded memory, meaning virtually
#  any number of arrays can be analyzed on also very limited computer
#  systems.
# }
#
# \usage{
#   @usage doACNE,AffymetrixCelSet
#   @usage doACNE,default
# }
#
# \arguments{
#  \item{csR, dataSet}{An @see "AffymetrixCelSet" (or the name of an
#   @see "AffymetrixCelSet").}
#  \item{fln}{If @TRUE, CRMAv2-style PCR fragment-length normalization
#   is performed, otherwise not.}
#  \item{drop}{If @TRUE, the RMA summaries are returned, otherwise
#   a named @list of all intermediate and final results.}
#  \item{verbose}{See @see "Verbose".}
#  \item{...}{Additional arguments used to set up @see "AffymetrixCelSet"
#   (when argument \code{dataSet} is specified).}
# }
#
# \value{
#   Returns a named @list, iff \code{drop == FALSE}, otherwise
#   a named @list of @see "aroma.core::AromaUnitTotalCnBinarySet"
#   and @see "aroma.core::AromaUnitFracBCnBinarySet".
# }
#
# \references{
#  [1] @include "../incl/OrtizM_etal_2010.Rd" \cr
# }
#
# @author "HB"
#*/###########################################################################
setMethodS3("doACNE", "AffymetrixCelSet", function(csR, fln=FALSE, drop=TRUE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csR':
  className <- "AffymetrixCelSet";
  if (!inherits(csR, className)) {
    throw(sprintf("Argument 'csR' is not a %s: %s", className, class(csR)[1]));
  }

  # Argument 'fln':
  fln <- Arguments$getVerbose(fln);

  # Argument 'drop':
  drop <- Arguments$getLogical(drop);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "ACNE");

  verbose && cat(verbose, "Data set");
  verbose && print(verbose, csR);

  # List of objects to be returned
  res <- list();
  if (!drop) {
    res <- c(res, list(csR=csR));
  }

  verbose && enter(verbose, "ACNE/CRMAv2/Allelic crosstalk calibration");
  acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2");
  verbose && print(verbose, acc);
  csC <- process(acc, verbose=verbose);
  verbose && print(verbose, csC);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(acc=acc, csC=csC));
  }

  # Clean up
  csR <- acc <- NULL;

  verbose && enter(verbose, "ACNE/CRMAv2/Base position normalization");
  bpn <- BasePositionNormalization(csC, target="zero");
  verbose && print(verbose, bpn);
  csN <- process(bpn, verbose=verbose);
  verbose && print(verbose, csN);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(bpn=bpn, csN=csN));
  }

  # Clean up
  csC <- bpn <- NULL;

  verbose && enter(verbose, "ACNE/Probe summarization");
  plm <- NmfSnpPlm(csN, mergeStrands=TRUE);
  verbose && print(verbose, plm);
  if (length(findUnitsTodo(plm)) > 0L) {
    # Fit CN probes quickly (~5-10s/array + some overhead)
    units <- fitCnProbes(plm, verbose=verbose);
    verbose && str(verbose, units);
    # Fit remaining units, i.e. SNPs (~5-10min/array)
    units <- fit(plm, verbose=verbose);
    verbose && str(verbose, units);
    units <- NULL;
  }
  # Clean up
  csN <- NULL;
  ces <- getChipEffectSet(plm);
  verbose && print(verbose, ces);
  verbose && exit(verbose);

  if (!drop) {
    res <- c(res, list(plm=plm));
  }

  # Clean up
  plm <- NULL;


  # PCR fragment-length normalization?
  if (fln) {
    verbose && enter(verbose, "ACNE/CRMAv2/PCR fragment-length normalization");
    fln <- FragmentLengthNormalization(ces, target="zero");
    verbose && print(verbose, fln);
    cesN <- process(fln, verbose=verbose);
    verbose && print(verbose, cesN);
    verbose && exit(verbose);

    if (!drop) {
      res <- c(res, list(fln=fln, cesN=cesN));
    }

    # Clean up
    fln <- ces <- NULL;
  } else {
    cesN <- ces;

    if (!drop) {
      res <- c(res, list(cesN=cesN));
    }
  }

  verbose && enter(verbose, "ACNE/Export to technology-independent data files");
  dsNList <- exportTotalAndFracB(cesN, verbose=verbose);
  verbose && print(verbose, dsNList);
  verbose && exit(verbose);

  # Clean up
  cesN <- NULL;

  if (!drop) {
    res <- c(res, list(dsNList=dsNList));
  }

  # Return only the final results?
  if (drop) {
    res <- dsNList;
  }

  verbose && exit(verbose);

  res;
}) # doACNE()


setMethodS3("doACNE", "default", function(dataSet, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "ACNE");

  verbose && enter(verbose, "ACNE/Setting up CEL set");
  csR <- AffymetrixCelSet$byName(dataSet, ..., verbose=less(verbose, 50),
                                                  .onUnknownArgs="ignore");
  verbose && print(verbose, csR);
  verbose && exit(verbose);

  res <- doACNE(csR, ..., verbose=verbose);

  # Clean up
  csR <- NULL;

  verbose && exit(verbose);

  res;
}) # doACNE()


############################################################################
# HISTORY:
# 2013-10-17
# o CLEANUP: Removed all explicit calls to gc().
# o CLEANUP: Dropped argument 'ram' to fit() of doACNE().
# o Turned doACNE() for character into a default method.
# o Created from doCRMAv1() in aroma.affymetrix.
############################################################################
