###########################################################################/**
# @RdocClass AromaUnitSignalBinarySet
#
# @title "The AromaUnitSignalBinarySet class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitSignalBinarySet object represents a set of
#  @see "AromaUnitSignalBinaryFile"s with \emph{identical} chip types.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaTabularBinarySet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("AromaUnitSignalBinarySet", function(...) {
  extend(AromaTabularBinarySet(...), "AromaUnitSignalBinarySet");
})


setMethodS3("findByName", "AromaUnitSignalBinarySet", function(static, ..., chipType=NULL) {
  NextMethod("findByName", subdirs=chipType);
}, static=TRUE, protected=TRUE)


setMethodS3("byName", "AromaUnitSignalBinarySet", function(static, name, tags=NULL, ..., chipType=NULL, paths=NULL, pattern="[.]asb$") {
  suppressWarnings({
    path <- findByName(static, name=name, tags=tags, chipType=chipType,
                                           ..., paths=paths, mustExist=TRUE);
  })

  suppressWarnings({
    byPath(static, path=path, ..., pattern=pattern);
  })
}, static=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN Interface API?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("validate", "AromaUnitSignalBinarySet", function(this, ...) {
  chipTypes <- lapply(this, FUN=getChipType);
  chipTypes <- unique(chipTypes);
  if (length(chipTypes) > 1) {
    throw("The located ", class(this)[1], " contains files with different chip types: ", paste(chipTypes, collapse=", "));
  }

  NextMethod("validate");
}, protected=TRUE)


setMethodS3("getPlatform", "AromaUnitSignalBinarySet", function(this, ...) {
  getPlatform(getOneFile(this), ...);
})


setMethodS3("getChipType", "AromaUnitSignalBinarySet", function(this, ...) {
  getChipType(getOneFile(this), ...);
})

setMethodS3("getAromaUgpFile", "AromaUnitSignalBinarySet", function(this, ...) {
  getAromaUgpFile(getOneFile(this), ...);
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END Interface API?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


############################################################################
# HISTORY:
# 2009-01-05
# o Renamed from AromaSignalBinarySet to AromaUnitSignalBinarySet.
# 2008-05-11
# o Created.
############################################################################
