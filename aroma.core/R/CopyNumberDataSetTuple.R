setMethodS3("getAverageFile", "AromaMicroarrayDataSet", abstract=TRUE);

# AD HOC. This is currently only implemented in the AffymetrixCelFile class,
# but should be implemented by others too, alternatively be replaced by
# a "better" method. /HB 2009-11-18.
setMethodS3("getXAM", "AromaMicroarrayDataFile", protected=TRUE, abstract=TRUE);


setConstructorS3("CopyNumberDataFile", function(...) {
  extend(Interface(), "CopyNumberDataFile");
})

setMethodS3("as.CopyNumberDataFile", "CopyNumberDataFile", function(this, ...) {
  this;
})


setMethodS3("hasAlleleBFractions", "CopyNumberDataFile", abstract=TRUE);
setMethodS3("hasStrandiness", "CopyNumberDataFile", abstract=TRUE);


setMethodS3("getNumberOfFilesAveraged", "CopyNumberDataFile", protected=TRUE, abstract=TRUE);


setConstructorS3("CopyNumberDataSet", function(...) {
  extend(Interface(), "CopyNumberDataSet");
})

setMethodS3("as.CopyNumberDataSet", "CopyNumberDataSet", function(this, ...) {
  this;
})

setMethodS3("hasAlleleBFractions", "CopyNumberDataSet", function(this, ...) {
  if (length(this) == 0L) {
    throw("Cannot infer hasAlleleBFractions(). No data files: ", getFullName(this));
  }
  hasAlleleBFractions(getOneFile(this));
})

setMethodS3("hasStrandiness", "CopyNumberDataSet", function(this, ...) {
  if (length(this) == 0L) {
    throw("Cannot infer hasStrandiness(). No data files: ", getFullName(this));
  }
  hasStrandiness(getOneFile(this));
})



setConstructorS3("CopyNumberDataSetTuple", function(...) {
  extend(Interface(), "CopyNumberDataSetTuple");
})

setMethodS3("as.CopyNumberDataSetTuple", "CopyNumberDataSetTuple", function(this, ...) {
  this;
})

setMethodS3("as.CopyNumberDataSetTuple", "list", function(this, ...) {
  n <- length(this);
  if (n == 0) {
    throw("Cannot coerce to CopyNumberDataSetTuple: The list is empty.");
  }

  # Infer the data set tuple class by coercing the first data set.
  set <- this[[1]];
  tuple <- as.CopyNumberDataSetTuple(set);

  # Then coerce the list  
  newInstance(tuple, this);
})

setMethodS3("hasAlleleBFractions", "CopyNumberDataSetTuple", function(this, ...) {
  cesList <- getSets(this);
  res <- sapply(cesList, FUN=hasAlleleBFractions);

  # Sanity check
  if (length(unique(res)) != 1) {
    throw("Inconsistent data sets: The ", class(this)[1], " contains data sets where some have combined the alleles and some have not.");
  }

  res <- res[1];
  res;
})

setMethodS3("hasStrandiness", "CopyNumberDataSetTuple", function(this, ...) {
  cesList <- getSets(this);
  res <- sapply(cesList, FUN=hasStrandiness);

  # Sanity check
  if (length(unique(res)) != 1) {
    throw("Inconsistent data sets: The ", class(this)[1], " contains data sets where some have data by strand and some have not.");
  }

  res <- res[1];
  res;
}) 



##############################################################################
# HISTORY:
# 2009-11-16
# o Created.
##############################################################################
