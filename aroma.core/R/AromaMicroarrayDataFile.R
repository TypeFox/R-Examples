###########################################################################/**
# @RdocClass AromaMicroarrayDataFile
#
# @title "The abstract AromaMicroarrayDataFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaMicroarrayDataFile object represents a single microarray data
#  file. Each such file originates from a specific chip type on a given
#  platform.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#
# \seealso{
#   An object of this class is typically part of an 
#   @see "AromaMicroarrayDataSet".
# }
#*/###########################################################################
setConstructorS3("AromaMicroarrayDataFile", function(...) {
  extend(GenericDataFile(...), c("AromaMicroarrayDataFile", 
                                               uses("FileCacheKeyInterface"))
  );
}, abstract=TRUE)


setMethodS3("getPlatform", "AromaMicroarrayDataFile", abstract=TRUE);


setMethodS3("getChipType", "AromaMicroarrayDataFile", abstract=TRUE);

 
setMethodS3("isAverageFile", "AromaMicroarrayDataFile", function(this, ...) {
  name <- getName(this);
  res <- (regexpr("^[.]average-", name) != -1);
  res;
})


############################################################################
# HISTORY:
# 2012-11-16
# o CLEANUP: Dropped (get|set)Label() for AromaMicroarrayDataFile.
# 2009-11-18
# o Added isAverageFile() for AromaMicroarrayDataFile.
# 2007-09-16
# o Created from AffymetrixFile.R.
############################################################################
