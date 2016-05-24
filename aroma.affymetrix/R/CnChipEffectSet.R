###########################################################################/**
# @RdocClass CnChipEffectSet
#
# @title "The CnChipEffectSet class"
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
#   \item{...}{Arguments passed to @see "SnpChipEffectSet".}
#   \item{combineAlleles}{A @logical indicating if the signals from
#      allele A and allele B are combined or not.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("CnChipEffectSet", function(..., combineAlleles="byFirstFile") {
  this <- extend(SnpChipEffectSet(...), c("CnChipEffectSet", uses("CopyNumberDataSet")));
  setCombineAlleles(this, combineAlleles);
  this;
})

setMethodS3("hasAlleleBFractions", "CnChipEffectSet", function(this, ...) {
  res <- (!this$combineAlleles);
  res;
})

setMethodS3("hasStrandiness", "CnChipEffectSet", function(this, ...) {
  res <- (!this$mergeStrands);
  res;
})


setMethodS3("byPath", "CnChipEffectSet", function(static, ..., combineAlleles="auto") {
  NextMethod("byPath", combineAlleles=combineAlleles);
}, static=TRUE, protected=TRUE)


setMethodS3("getAverageFile", "CnChipEffectSet", function(this, ...) {
  res <- NextMethod("getAverageFile");
  res$combineAlleles <- getCombineAlleles(this);
  res;
})

setMethodS3("getChipEffectFileClass", "CnChipEffectSet", function(static, ...) {
  CnChipEffectFile;
}, static=TRUE, private=TRUE)


setMethodS3("getCombineAlleles", "CnChipEffectSet", function(this, ...) {
  if (length(this) == 0L)
    return(FALSE);
  ce <- getOneFile(this);
  ce$combineAlleles;
})




setMethodS3("setCombineAlleles", "CnChipEffectSet", function(this, status, ...) {
  if (length(this) == 0L)
    return(FALSE);

  ce <- getOneFile(this);
  oldStatus <- ce$combineAlleles;

  if (identical(status, "byFirstFile")) {
    status <- oldStatus;
  }

  if (identical(status, "auto")) {
    status <- inferParameters(this, ...)$combineAlleles;
  }

  status <- Arguments$getLogical(status);
  lapply(this, FUN=function(ce) {
    ce$combineAlleles <- status;
  })
  invisible(oldStatus);
})


setMethodS3("inferParameters", "CnChipEffectSet", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Infer (mergeStrands, combineAlleles) parameters from stored data in quartet units in CEL set");

  # Identify units with quartets
  cdf <- getCdf(this);
  cdfPathname <- getPathname(cdf);
  nbrOfUnits <- .readCdfHeader(cdfPathname)$nunits;
  allUnits <- seq_len(nbrOfUnits);

  ce <- getOneFile(this);
  cePathname <- getPathname(ce);

  verbose && cat(verbose, "Pathname: ", cePathname);

  mergeStrands <- combineAlleles <- NA;
  while(length(allUnits) > 0) {
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
    units <- units[unitSizes == 4];

    if (length(units) > 0) {
      verbose && cat(verbose, "Scanning units:");
      verbose && str(verbose, units);
      # Infer parameters from 'intensities'
      values <- .readCelUnits(cePathname, units=units,
               readIntensities=TRUE, readStdvs=FALSE, dropArrayDim=TRUE);
      # Put quartets by columns
      values <- matrix(unlist(values, use.names=FALSE), nrow=4);
      # Keep only estimated units
      csums <- colSums(values);
      values <- values[,is.finite(csums) & (csums > 0),drop=FALSE];
      # Not needed anymore
      csums <- NULL;
      verbose && cat(verbose, "Values quartets:");
      verbose && print(verbose, values[,seq_len(min(ncol(values),6)),drop=FALSE]);
      if (ncol(values) > 0) {
        t <- .rowMeans(values);
        if (length(t) > 0) {
          isZero <- isZero(t);
          if (!all(isZero)) {
            combineAlleles <- isZero[2] && isZero[4];
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

  if (is.na(mergeStrands) || is.na(combineAlleles)) {
    throw("Failed to infer parameters 'mergeStrands' and 'combineAlleles' from chip-effect file: ", cePathname);
  }

  res <- list(combineAlleles=combineAlleles, mergeStrands=mergeStrands);

  verbose && str(verbose, res);
  verbose && exit(verbose);

  res;
}, private=TRUE)


setMethodS3("as.CopyNumberDataSetTuple", "CnChipEffectSet", function(this, ...) {
  CnChipEffectSetTuple(this, ...);
})


############################################################################
# HISTORY:
# 2009-11-18
# o Added as.CopyNumberDataSetTuple().
# 2009-11-17
# o Now CnChipEffectSet implements CopyNumberDataSet.
# 2008-05-16
# o Added support for setCombineAlleles(..., "byFirstFile").
# 2008-05-08
# o Made fromFiles() protected.
# 2007-11-20
# o BUG FIX: inferParams() would load all units if no units of size four
#   was found, because units <- units[unitsSizes == 4] => units == NULL.
# 2007-09-12
# o BUG FIX: A typo in setCombineAlleles() cause the function to always
#   return NULL instead of the last value.
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
# o Added fromFiles() which now infers 'combineAlleles' from the files.
# o Added inferParameters().
# 2006-11-22
# o Now getAverageFile() finally sets 'combineAlleles'.
# 2006-09-11
# o Created.
############################################################################
