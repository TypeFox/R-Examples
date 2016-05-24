# @author "HB"
setConstructorS3("AffymetrixPlatform", function(...) {
  extend(AromaPlatform(), "AffymetrixPlatform");
})

setMethodS3("getName", "AffymetrixPlatform", function(static, ...) {
  "Affymetrix";
})

setMethodS3("getUnitNamesFile", "AffymetrixPlatform", function(static, ...) {
  AffymetrixCdfFile$byName(...);
}, static=TRUE)

setMethodS3("getUnitTypesFile", "AffymetrixPlatform", function(static, ...) {
  AffymetrixCdfFile$byName(...);
}, static=TRUE)




############################################################################
# HISTORY:
# 2009-07-08
# o Added getUnitTypesFile() for AffymetrixPlatform.
# 2008-05-18
# o Created.
############################################################################
