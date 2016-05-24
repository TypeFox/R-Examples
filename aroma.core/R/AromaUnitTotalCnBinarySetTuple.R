setConstructorS3("AromaUnitTotalCnBinarySetTuple", function(..., .setClass="AromaUnitTotalCnBinarySet") {
  extend(AromaMicroarrayDataSetTuple(..., .setClass=.setClass), c("AromaUnitTotalCnBinarySetTuple", uses("CopyNumberDataSetTuple")));
})


setMethodS3("as.AromaUnitTotalCnBinarySetTuple", "AromaUnitTotalCnBinarySetTuple", function(this, ...) {
  # Nothing do to
  this;
})


setMethodS3("getListOfUnitNamesFiles", "AromaUnitTotalCnBinarySetTuple", function(this, ...) {
  sets <- getSets(this);
  lapply(sets, FUN=getUnitNamesFile);
})


setMethodS3("getUnitNamesFile", "AromaUnitTotalCnBinarySet", function(this, ...) {
  if (length(this) == 0) {
    throw("Cannot locate unit names file. Data set is empty: ", getFullName(this));
  }

  getUnitNamesFile(getOneFile(this), ...);
})


###########################################################################
# HISTORY:
# 2009-11-19
# o Created.
###########################################################################
