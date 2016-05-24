###########################################################################/**
# @RdocClass AromaUnitPscnBinaryFile
#
# @title "The AromaUnitPscnBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitPscnBinaryFile is a @see "AromaUnitSignalBinaryFile"
#  that holds total copy number signals (TCNs) and allele B fractions (BAFs).
#  The TCNs can either be on an unknown scale or ratios relative
#  to a reference.  The signals are always stored on the original scale,
#  i.e. they are never stored on the logaritmic scale.
#  The BAFs are always on a [0-eps,1+eps] scale, where eps >= 0.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaUnitSignalBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/########################################################################### 
setConstructorS3("AromaUnitPscnBinaryFile", function(...) {
  extend(AromaUnitSignalBinaryFile(...), c("AromaUnitPscnBinaryFile", uses("CopyNumberDataFile"))
  );
})


setMethodS3("hasTotalCopyNumberRatios", "AromaUnitPscnBinaryFile", function(this, ...) {
  hasTag(this, "ratio");
})


setMethodS3("hasAlleleBFractions", "AromaUnitPscnBinaryFile", function(this, ...) {
  # By definition, always TRUE
  TRUE;
})

setMethodS3("hasStrandiness", "AromaUnitPscnBinaryFile", function(this, ...) {
  # For now always FALSE, due to how the (super)classes are defined.
  # /HB 2009-11-19
  FALSE;
})
 
setMethodS3("getDefaultColumnNames", "AromaUnitPscnBinaryFile", function(this, ...) {
 stopifnot(nbrOfColumns(this) == 2L);
 c("total", "fracB");
}, protected=TRUE)


setMethodS3("extractRawCopyNumbers", "AromaUnitPscnBinaryFile", function(this, ..., logBase=2L, clazz=RawCopyNumbers) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'logBase':
  if (!is.null(logBase)) {
    logBase <- Arguments$getDouble(logBase, range=c(1, 10));
  }

  # Total CN signals are always stored on the original scale in column #1.
  cn <- extractRawGenomicSignals(this, column=1L, ..., clazz=clazz);

  # Convert to the correct logarithmic base
  cn <- extractRawCopyNumbers(cn, logBase=logBase);

  cn;
})


setMethodS3("allocate", "AromaUnitPscnBinaryFile", function(static, ..., platform, chipType, types=c("double", "double"), sizes=c(4L, 4L), signed=c(TRUE, TRUE), footer=list()) {
  # Argument 'platform':
  platform <- Arguments$getCharacter(platform, length=c(1,1));

  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  # Create tabular binary file
  res <- NextMethod("allocate", generic="allocate", types=types, sizes=sizes, signeds=signed);


  # Write attributes to footer
  attrs <- list(
    createdOn=format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE),
    platform=platform, 
    chipType=chipType
  );
  footer <- c(attrs, footer);
  writeFooter(res, footer);

  res;
}, static=TRUE, protected=TRUE)


############################################################################
# HISTORY:
# 2012-07-21
# o Added allocate().
# 2012-07-20
# o Created from AromaUnitTotalCnBinaryFile.R.
############################################################################
