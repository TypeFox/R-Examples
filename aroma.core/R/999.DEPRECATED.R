## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## DEFUNCT
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 2016-01-05: Defunct
## 2015-09-18: Deprecated
setMethodS3("apply", "SampleAnnotationFile", function(...) {
  .Defunct("applyTo")
  applyTo(...)
}, protected=TRUE)


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## DEPRECATED
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("whatDataType", "default", function(type, ...) {
  .Deprecated(msg="whatDataType() is deprecated (as of aroma.core 2.14.0) with no alternative implementation. Please contact the maintainer of the aroma.core package if you wish that this function should remain available.")

  if (type == "byte") {
    what <- "integer"
    size <- 1
  } else if (type == "short") {
    what <- "integer"
    size <- 2
  } else if (type == "integer") {
    what <- "integer"
    size <- 4
  } else if (type == "long") {
    what <- "integer"
    size <- 8
  } else if (type == "float") {
    what <- "double"
    size <- 4
  } else if (type == "double") {
    what <- "double"
    size <- 8
  } else if (type == "logical") {
    what <- "integer"
    size <- 1
  } else if (type == "raw") {
    what <- "raw"
    size <- 1
  } else {
    what <- NA
    size <- NA
  }

  list(what=what, size=size)
}, private=TRUE, deprecated=TRUE)



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## TO DEPRECATE
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



############################################################################
# HISTORY:
# 2014-02-28
# o CLEANUP: Defuncted previously deprecated downloadPackagePatch() and
#   patchPackage() as well as patch() for AromaPackage.
# 2014-02-28
# o CLEANUP: Removed defunct methods, i.e. is(Homo|Hetero)zygote() and
#   getPhysicalPositions().
# 2013-04-20
# o CLEANUP: Made is(Homo|Hetero)zygote() & getPhysicalPositions() defunct.
# 2012-10-14
# o Created 999.DEPRECATED.R.
# 2011-02-18
# o CLEANUP: Deprecated static method importFromTable() for FileMatrix.
############################################################################
