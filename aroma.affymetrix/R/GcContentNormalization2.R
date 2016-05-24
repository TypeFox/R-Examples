###########################################################################/**
# @RdocClass GcContentNormalization2
#
# @title "The GcContentNormalization2 class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization method that corrects for
#  annotation-data covariate effects on copy-number chip-effect estimates.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "AdditiveCovariatesNormalization".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("GcContentNormalization2", function(...) {
  extend(AdditiveCovariatesNormalization(...), "GcContentNormalization2");
})


setMethodS3("getCovariates", "GcContentNormalization2", function(this, units=NULL, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'units':
  cdf <- getCdf(this);
  if (!is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    units <- Arguments$getIndices(units, max=nbrOfUnits);
  }


  verbose && enter(verbose, "Retrieving GC content covariates");

  chipType <- getChipType(cdf);
  chipType <- gsub(",monocell", "", chipType);
  verbose && cat(verbose, "Chip type: ", chipType);
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);

  X <- NULL;
  # Try 1: Use an unit GC content (UGC) file
  tryCatch({
    ugc <- AromaUgcFile$byChipType(chipType);
    X <- ugc[units,1,drop=TRUE];
    X <- as.matrix(X);
  }, error = function(ex) {
  })

  # Try 2: Use a TSV file (deprecated; kept for backward compatibility)
  if (is.null(X)) {
    tryCatch({
      chipTypeS <- gsub(",.*", "", chipType);
      tsv <- AffymetrixTsvFile$byChipType(chipTypeS);
      X <- getGc(tsv, units=units);
      X <- as.matrix(X);
    }, error = function(ex) {
    })
  }

  if (is.null(X)) {
    throw("Failed to retrieve covariates. No GC-content annotation file found: ", chipType);
  }

  verbose && cat(verbose, "GC contents:");
  verbose && str(verbose, X);

  verbose && exit(verbose);

  X;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2009-03-22
# o Created.
############################################################################
