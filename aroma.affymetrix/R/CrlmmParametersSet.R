###########################################################################/**
# @RdocClass CrlmmParametersSet
#
# @title "The CrlmmParametersSet class"
#
# \description{
#  @classhierarchy
#
#  An CrlmmParametersSet object represents a set of
#  @see "CrlmmParametersFile"s with \emph{identical} chip types.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to
#     @see "aroma.core::AromaUnitSignalBinarySet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("CrlmmParametersSet", function(...) {
  extend(AromaUnitSignalBinarySet(...), "CrlmmParametersSet");
})


setMethodS3("byName", "CrlmmParametersSet", function(static, name, tags=NULL, ..., chipType=NULL, paths="crlmmData(|,.*)/") {
  suppressWarnings({
    path <- findByName(static, name=name, tags=tags, chipType=chipType,
                                           ..., paths=paths, mustExist=TRUE);
  })

  byPath(static, path=path, ...);
}, static=TRUE)

setMethodS3("byPath", "CrlmmParametersSet", function(static, ...) {
  suppressWarnings({
    NextMethod("byPath", pattern=".*,CRLMM[.]atb$$");
  })
})


setMethodS3("findUnitsTodo", "CrlmmParametersSet", function(this, ...) {
  # Look into the chip-effect file that comes last in a lexicographic
  # order, becuase that is updated last.
  names <- getFullNames(this);
  idx <- order(names, decreasing=TRUE)[1];
  df <- this[[idx]];
  findUnitsTodo(df, ...);
})



############################################################################
# HISTORY:
# 2011-02-24
# o Expanded the searched root paths to be crlmmData(|,.*)/.
# 2010-05-08
# o Now all findUnitsTodo() for data sets checks the data file that comes
#   last in a lexicographic ordering.  This is now consistent with how
#   the summarization methods updates the files.  Before it was use to be
#   the one that is last in the data set.
# 2008-12-08
# o Added findUnitsTodo() and extractCalls().
# 2008-12-06
# o Created.
############################################################################
