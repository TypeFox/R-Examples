###########################################################################/**
# @RdocClass SnpChipEffectSet
#
# @title "The SnpChipEffectSet class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in the probe-level models.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ChipEffectSet".}
#   \item{mergeStrands}{Specifies if the strands are merged or not for these
#      estimates.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("SnpChipEffectSet", function(..., mergeStrands="byFirstFile") {
  this <- extend(ChipEffectSet(...), "SnpChipEffectSet");
  setMergeStrands(this, mergeStrands);
  this;
})


setMethodS3("byPath", "SnpChipEffectSet", function(static, ..., mergeStrands="auto") {
  NextMethod("byPath", mergeStrands=mergeStrands);
}, static=TRUE, protected=TRUE)



setMethodS3("getAverageFile", "SnpChipEffectSet", function(this, ...) {
  res <- NextMethod("getAverageFile");
  res$mergeStrands <- getMergeStrands(this);
  res;
})



setMethodS3("getChipEffectFileClass", "SnpChipEffectSet", function(static, ...) {
  SnpChipEffectFile;
}, static=TRUE, private=TRUE)

setMethodS3("getMergeStrands", "SnpChipEffectSet", function(this, ...) {
  if (length(this) == 0L)
    return(FALSE);
  ce <- getOneFile(this);
  ce$mergeStrands;
})

setMethodS3("setMergeStrands", "SnpChipEffectSet", function(this, status, ...) {
  if (length(this) == 0)
    return(FALSE);

  oldStatus <- getMergeStrands(this);

  if (identical(status, "byFirstFile")) {
    status <- oldStatus;
  }

  if (identical(status, "auto")) {
    status <- inferParameters(this, ...)$mergeStrands;
  }

  # Argument 'status':
  status <- Arguments$getLogical(status);

  # Update all chip-effect files
  lapply(this, FUN=function(ce) {
    ce$mergeStrands <- status;
  })

  invisible(oldStatus);
})


setMethodS3("inferParameters", "SnpChipEffectSet", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Infer (mergeStrands) parameters from stored data in quartet units in CEL set");

  # Identify units with quartets
  cdf <- getCdf(this);
  cdfPathname <- getPathname(cdf);
  nbrOfUnits <- .readCdfHeader(cdfPathname)$nunits;
  allUnits <- seq_len(nbrOfUnits);

  ce <- getOneFile(this, mustExist=TRUE);
  cePathname <- getPathname(ce);

  verbose && cat(verbose, "Pathname: ", cePathname);

  mergeStrands <- NA;
  while(length(allUnits) > 0L) {
    verbose && cat(verbose, "Units left: ", length(allUnits));
    uu <- seq_len(min(10e3,length(allUnits)));
    units <- allUnits[uu];
    allUnits <- allUnits[-uu];
    # Not needed anymore
    uu <- NULL;

    # Identify units that are quartets
    unitSizes <- .readCdfGroupNames(cdfPathname, units=units);
    names(unitSizes) <- NULL;
    unitSizes <- sapply(unitSizes, FUN=length);
    units <- units[unitSizes == 4L];
    # Not needed anymore
    unitSizes <- NULL;
    if (length(units) > 0L) {
      verbose && cat(verbose, "Scanning units:");
      verbose && str(verbose, units);
      # Infer parameters from 'intensities'
      values <- .readCelUnits(cePathname, units=units,
                   readIntensities=TRUE, readStdvs=FALSE, dropArrayDim=TRUE);
      # Not needed anymore
      units <- NULL;

      # Put quartets by columns
      values <- matrix(unlist(values, use.names=FALSE), nrow=4L);

      # Keep only estimated units without NAs
      csums <- colSums(values);
      values <- values[,is.finite(csums) & (csums > 0),drop=FALSE];
      # Not needed anymore
      csums <- NULL;
      verbose && cat(verbose, "Values quartets:");
      verbose && print(verbose, values[,seq_len(min(ncol(values),6)),drop=FALSE]);
      if (ncol(values) > 0L) {
        t <- .rowMeans(values);
        if (length(t) > 0L) {
          isZero <- isZero(t);
          if (!all(isZero)) {
            mergeStrands <- isZero[3] && isZero[4];
            break;
          }
          # Not needed anymore
          isZero <- NULL;
        }
        # Not needed anymore
        t <- NULL;
      }
      # Not needed anymore
      values <- NULL;
    }
  } # while(...)

  if (is.na(mergeStrands)) {
    throw("Failed to infer parameter 'mergeStrands' from chip-effect file: ", cePathname);
  }

  res <- list(mergeStrands=mergeStrands);

  verbose && str(verbose, res);
  verbose && exit(verbose);

  res;
}, private=TRUE)



############################################################################
# HISTORY:
# 2008-05-16
# o Added support for setMergeStrands(..., "byFirstFile").
# 2008-05-08
# o Made fromFiles() protected.
# 2007-11-20
# o BUG FIX: inferParams() would load all units if no units of size four
#   was found, because units <- units[unitsSizes == 4] => units == NULL.
# 2007-03-23
# o Now inferParameters() are looking at the 'intensity' (==theta) field
#   instead of 'stdvs'.  The reason for this is that 'stdvs' might be all
#   zeros, e.g. after a fragment-length normalization.
# 2007-02-20
# o BUG FIX: inferParameters() would give an error if some estimates were
#   NAs.
# 2007-02-19
# o BUG FIX: If inferParameters() where called on a chip-effect set where
#   some units where not yet estimated, an error would be generated.
# 2007-01-11
# o Added fromFiles() which now infers 'mergeStrands' from the files.
# o Added inferParameters().
# 2006-11-22
# o Now getAverageFile() finally sets 'mergeStrands'.
# 2006-09-11
# o Created.
############################################################################
