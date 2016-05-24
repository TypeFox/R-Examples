###########################################################################/**
# @RdocClass QualityAssessmentFile
#
# @title "The QualityAssessmentFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents probe-level QC information (residuals, weights, etc.)
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixCelFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "KS"
#
# \seealso{
#   An object of this class is typically part of a @see "QualityAssessmentSet".
# }
#*/###########################################################################
setConstructorS3("QualityAssessmentFile", function(...) {
  this <- extend(AffymetrixCelFile(...), "QualityAssessmentFile",
    "cached:.firstCells" = NULL
  );

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("findUnitsTodo", "QualityAssessmentFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Identifying non-assigned units in QC file");

  verbose && cat(verbose, "Pathname: ", getPathname(this));
  if (is.null(units)) {
    units <- seq_len(nbrOfUnits(getCdf(this)));
  }

  # Read 'pixels' from each unit
  verbose && enter(verbose, "Reading data for these ", length(units), " units");
  value <- .readCelUnits(getPathname(this), units=units, readIntensities=FALSE,
                        readStdvs=FALSE, readPixels=TRUE, dropArrayDim=TRUE);
  verbose && exit(verbose);

  verbose && enter(verbose, "Looking for pixels == 0 indicating non-assigned units");
  # Identify units for which all pixels == 0.
  allZeroPixels <- sapply(value, function(x) {
    all(x[[1]][[1]]==0)
  }, USE.NAMES=FALSE);
  value <- which(allZeroPixels);
  if (!is.null(units))
    value <- units[value];
  verbose && str(verbose, value);
  verbose && exit(verbose);

  verbose && exit(verbose);

  value;
})


############################################################################
# HISTORY:
# 2007-02-12 /HB
# o Updated the verbose output for findUnitsTodo().
# 2007-01-12
# o Created.
############################################################################
