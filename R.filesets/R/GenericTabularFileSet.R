###########################################################################/**
# @RdocClass GenericTabularFileSet
#
# @title "The GenericTabularFileSet class"
#
# \description{
#  @classhierarchy
#
#  An GenericTabularFileSet object represents a set of
#  @see "GenericTabularFile"s.
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
setConstructorS3("GenericTabularFileSet", function(...) {
  extend(GenericDataFileSet(...), "GenericTabularFileSet");
})



setMethodS3("extractMatrix", "GenericTabularFileSet", function(this, files=NULL, ..., drop=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfFiles <- length(this);

  # Argument 'files': [DEPRECATED]
  if (is.null(files)) {
    files <- seq_len(nbrOfFiles);
  } else {
    .Deprecated("Argument 'files' of extractMatrix() for GenericTabularFileSet is deprecated. Use extractMatrix(ds[files], ...) instead.")
    files <- Arguments$getIndices(files, max=nbrOfFiles);
    nbrOfFiles <- length(files);
  }


  data <- NULL;
  for (kk in seq_len(nbrOfFiles)) {
    file <- files[kk];
    dataFile <- this[[file]];
    dataKK <- extractMatrix(dataFile, ...);

    if (is.null(data)) {
      naValue <- vector(storage.mode(dataKK), length=1);
      data <- matrix(naValue, nrow=nrow(dataKK), ncol=nbrOfFiles);
      colnames(data) <- getNames(this)[files];
    }

    data[,kk] <- dataKK;
    # Not needed anymore
    dataKK <- NULL;
  }

  # Drop singelton dimensions?
  if (drop) {
    data <- drop(data);
  }

  data;
})


############################################################################
# HISTORY:
# 2008-05-12
# o Created.
############################################################################
