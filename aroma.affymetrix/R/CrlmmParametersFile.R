###########################################################################/**
# @RdocClass CrlmmParametersFile
#
# @title "The CrlmmParametersFile class"
#
# \description{
#  @classhierarchy
#
#  An CrlmmParametersFile is a @see "aroma.core::AromaUnitSignalBinaryFile".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "aroma.core::AromaUnitSignalBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("CrlmmParametersFile", function(...) {
  extend(AromaUnitSignalBinaryFile(...), "CrlmmParametersFile"
  );
})


setMethodS3("allocate", "CrlmmParametersFile", function(static, ..., nbrOfStrands=2, types=rep("double", times=1+3*nbrOfStrands), sizes=rep(4L, times=1+3*nbrOfStrands), signed=rep(TRUE, times=1+3*nbrOfStrands)) {
  NextMethod("allocate", types=types, sizes=sizes, signed=signed);
})



setMethodS3("findUnitsTodo", "CrlmmParametersFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Identifying non-fitted units in file");
  verbose && cat(verbose, "Pathname: ", getPathname(this));

  # Reading all calls
  values <- this[,1,drop=TRUE];

  units <- which(values == 0);
  verbose && exit(verbose);

  units;
})


setMethodS3("readParameter", "CrlmmParametersFile", function(this, name, mode="character", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'name':
  name <- Arguments$getCharacter(name, length=c(1,1));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading parameter (stored in file footer)");
  footer <- readFooter(this);
  key <- "parameters";
  params <- footer[[key]];
  res <- params[[name]];

  storage.mode(res) <- mode;
  verbose && exit(verbose);

  res;
})



setMethodS3("updateParameter", "CrlmmParametersFile", function(this, name, value, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'name':
  name <- Arguments$getCharacter(name, length=c(1,1));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Updating parameter (stored in file footer)");
  footer <- readFooter(this);
  verbose && cat(verbose, "File footer before:");
  verbose && str(verbose, footer);

  key <- "parameters";
  params <- footer[[key]];
  if (is.null(params)) {
    params <- list();
  }
  params[[name]] <- value;
  footer[[key]] <- params;

  verbose && cat(verbose, "Updated footer:");
  verbose && str(verbose, footer);

  res <- writeFooter(this, footer);
  verbose && exit(verbose);

  invisible(res);
})



############################################################################
# HISTORY:
# 2009-01-12
# o Added read-/updateParameter().
# 2008-12-08
# o Added findUnitsTodo() and extractCalls().
# 2008-12-05
# o Created.
############################################################################
