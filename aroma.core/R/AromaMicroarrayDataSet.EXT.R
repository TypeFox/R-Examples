setMethodS3("as.AromaMicroarrayDataSetTuple", "AromaMicroarrayDataSet", function(this, ...) {
  classNames <- class(this);
  classNames <- sprintf("%sTuple", classNames);
  methodNames <- sprintf("as.%s", classNames);

  # Try first as.Nnn() method available
  for (kk in seq_along(methodNames)) {
    methodName <- methodNames[kk];
    if (exists(methodName, mode="function")) {
      fcn <- get(methodName, mode="function");
      res <- fcn(this, ...);
      return(res);
    }
  } # for (kk ...)

  throw("Failed to coerce ", class(this)[1], ", to an AromaMicroarrayDataSetTuple: ", getFullName(this));
})


setMethodS3("as.AromaMicroarrayDataSetList", "AromaMicroarrayDataSet", function(this, ...) {
  classNames <- class(this);
  classNames <- sprintf("%sList", classNames);
  methodNames <- sprintf("as.%s", classNames);

  # Try first as.Nnn() method available
  for (kk in seq_along(methodNames)) {
    methodName <- methodNames[kk];
    if (exists(methodName, mode="function")) {
      fcn <- get(methodName, mode="function");
      res <- fcn(this, ...);
      return(res);
    }
  } # for (kk ...)

  throw("Failed to coerce ", class(this)[1], ", to an AromaMicroarrayDataSetList: ", getFullName(this));
})


##############################################################################
# HISTORY:
# 2009-12-30
# o Added as.AromaMicroarrayDataSetList() to AromaMicroarrayDataSet.
# 2009-11-18
# o Created.
##############################################################################

