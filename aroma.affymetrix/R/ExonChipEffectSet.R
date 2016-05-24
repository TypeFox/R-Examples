###########################################################################/**
# @RdocClass ExonChipEffectSet
#
# @title "The ExonChipEffectSet class"
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
#   \item{mergeGroups}{Specifies if groups (individual exons in a CDF
#        file) are merged or not for these estimates, i.e. whether
#        transcript-level expression is to be estimated.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "KS, HB"
#*/###########################################################################
setConstructorS3("ExonChipEffectSet", function(..., mergeGroups=TRUE) {
  this <- extend(ChipEffectSet(...), "ExonChipEffectSet");
  setMergeGroups(this, mergeGroups);
  this;
})

setMethodS3("byPath", "ExonChipEffectSet", function(static, ..., mergeGroups="auto") {
  NextMethod("byPath", mergeGroups=mergeGroups);
}, static=TRUE, protected=TRUE)



setMethodS3("getAverageFile", "ExonChipEffectSet", function(this, ...) {
  res <- NextMethod("getAverageFile");
  res$mergeGroups <- getMergeGroups(this);
  res;
})



setMethodS3("getChipEffectFileClass", "ExonChipEffectSet", function(static, ...) {
  ExonChipEffectFile;
}, static=TRUE, private=TRUE)

setMethodS3("getMergeGroups", "ExonChipEffectSet", function(this, ...) {
  if (length(this) == 0L)
    return(FALSE);
  ce <- getOneFile(this);
  ce$mergeGroups;
})

setMethodS3("setMergeGroups", "ExonChipEffectSet", function(this, status, ...) {
  if (length(this) == 0)
    return(FALSE);

  oldStatus <- getMergeGroups(this);

#  if (identical(status, "auto"))
#    status <- inferParameters(this)$mergeGroups;

  # Argument 'status':
  status <- Arguments$getLogical(status);

  # Update all chip-effect files
  lapply(this, FUN=function(ce) {
    ce$mergeGroups <- status;
  })

  invisible(oldStatus);
})

setMethodS3("getFirstCellPerUnitIndices", "ExonChipEffectSet", function(this, ...) {

  cdf <- getCdf(this);
  idx <- getFirstCellIndices(cdf, ...);
  idx <- lapply(lapply(idx, FUN=.subset2, 1), FUN=.subset2, 1);
  idx <- unlist(idx, use.names=FALSE);
  idx;
}, protected=TRUE)


setMethodS3("findUnitsTodo", "ExonChipEffectSet", function(this, ...) {
  # Look into the chip-effect file that comes last in a lexicographic
  # order, becuase that is updated last.
  names <- getFullNames(this);
  idx <- order(names, decreasing=TRUE)[1];
  df <- this[[idx]];
  findUnitsTodo(df, ...);
})



############################################################################
# HISTORY:
# 2010-05-08
# o Now all findUnitsTodo() for data sets checks the data file that comes
#   last in a lexicographic ordering.  This is now consistent with how
#   the summarization methods updates the files.  Before it was use to be
#   the one that is last in the data set.
# 2008-05-08
# o Made fromFiles() protected.
# 2007-02-08
# o Created (based on SnpChipEffectSet.R following chat with HB on
#   2007-02-07).
############################################################################
