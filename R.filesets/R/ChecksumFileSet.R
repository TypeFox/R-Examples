###########################################################################/**
# @RdocClass ChecksumFileSet
#
# @title "The ChecksumFileSet class"
#
# \description{
#  @classhierarchy
#
#  An ChecksumFileSet object represents a set of @see "ChecksumFile"s.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("ChecksumFileSet", function(...) {
  extend(GenericDataFileSet(...), "ChecksumFileSet")
})


setMethodS3("getChecksumFileSet", "GenericDataFileSet", function(this, ...) {
  files <- vector("list", length=length(this));
  for (ii in seq_along(this)) {
    file <- this[[ii]];
    files[[ii]] <- getChecksumFile(file, ...);
  } # for (ii ...)
  ChecksumFileSet(files);
})


setMethodS3("validate", "ChecksumFileSet", function(this, ...) {
  lapply(this, FUN=validate, ...);
  invisible(TRUE);
})


setMethodS3("readChecksums", "ChecksumFileSet", function(ds, ...) {
  sapply(ds, FUN=readChecksum, ...)
})


############################################################################
# HISTORY:
# 2014-08-19
# o Added readChecksums() for ChecksumFileSet.
# 2013-11-19
# o Added ChecksumFile.
# o Created.
############################################################################
