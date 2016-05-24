###########################################################################/**
# @RdocClass AromaUnitFracBCnBinaryFile
#
# @title "The AromaUnitFracBCnBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitFracBCnBinaryFile is a @see "AromaUnitTabularBinaryFile".
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
setConstructorS3("AromaUnitFracBCnBinaryFile", function(...) {
  extend(AromaUnitSignalBinaryFile(...), "AromaUnitFracBCnBinaryFile"
  );
})


setMethodS3("extractRawAlleleBFractions", "AromaUnitFracBCnBinaryFile", function(this, ..., clazz=RawAlleleBFractions) {
  extractRawGenomicSignals(this, ..., clazz=clazz);
})



############################################################################
# HISTORY:
# 2009-05-17
# o Added extractRawAlleleBFractions() for AromaUnitFracBCnBinaryFile.
# 2008-05-11
# o Created.
############################################################################
