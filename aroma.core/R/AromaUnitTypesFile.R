###########################################################################/**
# @RdocClass AromaUnitTypesFile
#
# @title "The AromaUnitTypesFile class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaUnitTabularBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/########################################################################### 
setConstructorS3("AromaUnitTypesFile", function(...) {
  extend(AromaUnitSignalBinaryFile(...), c("AromaUnitTypesFile", 
                                                        uses("UnitTypesFile"))
  );
})

setMethodS3("getChipType", "AromaUnitTypesFile", function(this, ...) {
  readFooter(this)$chipType;
})

setMethodS3("getPlatform", "AromaUnitTypesFile", function(this, ...) {
  readFooter(this)$platform;
})


setMethodS3("getUnitTypes", "AromaUnitTypesFile", function(this, ...) {
  data <- extractMatrix(this, column=1, drop=TRUE, ...);
  ftr <- readFooter(this);
  types <- ftr$types;
  attr(data, "types") <- types;
  data;
})

setMethodS3("allocate", "AromaUnitTypesFile", function(static, ..., types=c("integer"), sizes=1L) { 
  NextMethod("allocate", types=types, sizes=sizes);
}, static=TRUE, protected=TRUE)


setMethodS3("importFromUnitTypesFile", "AromaUnitTypesFile", function(this, utf, ..., verbose=FALSE) { 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'utf':
  utf <- Arguments$getInstanceOf(utf, "UnitTypesFile");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Importing unit types from ", class(utf)[1]);
  verbose && print(verbose, utf);
  unitTypes <- getUnitTypes(utf, verbose=less(verbose, 10));
  this[,1] <- unitTypes;
  verbose && exit(verbose);

  invisible(TRUE);
}, static=TRUE)






############################################################################
# HISTORY:
# 2009-07-09
# o Created.
############################################################################
