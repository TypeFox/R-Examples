###########################################################################/**
# @RdocClass AromaUnitCallSet
#
# @title "The AromaUnitCallSet class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitCallSet object represents a set of
#  @see "AromaUnitCallFile"s with \emph{identical} chip types.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaUnitSignalBinarySet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("AromaUnitCallSet", function(...) {
  extend(AromaUnitSignalBinarySet(...), "AromaUnitCallSet");
})


setMethodS3("findByName", "AromaUnitCallSet", function(static, ..., paths="callData/") {
  NextMethod("findByName", paths=paths);
}, static=TRUE, protected=TRUE)


setMethodS3("byPath", "AromaUnitCallSet", function(static, ..., pattern=".*[.]acf$") {
  suppressWarnings({
    NextMethod("byPath", pattern=pattern);
  })
}, static=TRUE, protected=TRUE)


setMethodS3("findUnitsTodo", "AromaUnitCallSet", function(this, ...) {
  # Look into the last chip-effect file since that is updated last
  df <- this[[length(this)]];
  findUnitsTodo(df, ...);
})


setMethodS3("extractCallArray", "AromaUnitCallSet", function(this, ..., drop=FALSE) {
  res <- NULL;
  for (kk in seq_along(this)) {
    df <- this[[kk]];
    values <- extractCalls(df, ..., drop=FALSE);

    if (is.null(res)) {
      dim <- dim(values);
      dim[length(dim)] <- length(this);
      res <- array(values[1], dim=dim);
    }

    res[,,kk] <- values;
  } # for (kk ...)

  # Drop singletons?
  if (drop) {
    res <- drop(res);
  }

  res;
})


setMethodS3("extractCalls", "AromaUnitCallSet", function(this, ...) {
  extractCallArray(this, ...);
})


############################################################################
# HISTORY:
# 2009-01-04
# o Created.
############################################################################
