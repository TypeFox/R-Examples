###########################################################################/**
# @RdocClass AromaUnitTotalCnBinaryFile
#
# @title "The AromaUnitTotalCnBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitTotalCnBinaryFile is a @see "AromaUnitSignalBinaryFile".
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
setConstructorS3("AromaUnitTotalCnBinaryFile", function(...) {
  extend(AromaUnitSignalBinaryFile(...), c("AromaUnitTotalCnBinaryFile", uses("CopyNumberDataFile"))
  );
})


setMethodS3("hasAlleleBFractions", "AromaUnitTotalCnBinaryFile", function(this, ...) {
  # By definition, always FALSE
  FALSE;
})

setMethodS3("hasStrandiness", "AromaUnitTotalCnBinaryFile", function(this, ...) {
  # For now always FALSE, due to how the (super)classes are defined. 
  # /HB 2009-11-19
  FALSE;
})
 

setMethodS3("extractRawCopyNumbers", "AromaUnitTotalCnBinaryFile", function(this, ..., logBase=2, clazz=RawCopyNumbers) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'logBase':
  if (!is.null(logBase)) {
    logBase <- Arguments$getDouble(logBase, range=c(1, 10));
  }


  cn <- extractRawGenomicSignals(this, ..., clazz=clazz);

  logBase0 <- NULL;
  if (hasTag(this, "log2ratio")) {
    logBase0 <- 2;
  } else if (hasTag(this, "log10ratio")) {
    logBase0 <- 10;
  } else if (hasTag(this, "logRatio")) {
    # AD HOC. Now standard
    logBase0 <- 10;
  }

  cn <- setBasicField(cn, ".yLogBase", logBase0);

  # Convert to the correct logarithmic base
  cn <- extractRawCopyNumbers(cn, logBase=logBase);

  cn;
})



############################################################################
# HISTORY:
# 2012-07-20
# o CLEANUP: Moved getNumberOfFilesAveraged() to superclass.
# 2009-11-22
# o Added argument 'logBase=2' to extractRawCopyNumbers() for 
#   AromaUnitTotalCnBinaryFile.  Now it is possible to extract either
#   non-logged or logged CN ratios.
# 2009-11-20
# o Added getNumberOfFilesAveraged for of AromaUnitTotalCnBinaryFile.
# 2009-05-17
# o Now extractRawCopyNumbers() of AromaUnitTotalCnBinaryFile utilized
#   extractRawGenomicSignals() of the superclass.
# 2009-02-17
# o Added argument 'units' to extractRawCopyNumbers().
# 2009-02-16
# o Now extractRawCopyNumbers() also includes the full (sample) name.
# 2008-06-12
# o Now extractRawCopyNumbers() adds annotation data to the returned object,
#   i.e. platform, chipType, and fullname.
# 2008-05-21
# o Added extractRawCopyNumbers().
# 2008-05-11
# o Created.
############################################################################
