setConstructorS3("AromaUnitGcContentFile", function(...) {
  this <- extend(AromaUnitTabularBinaryFile(...), "AromaUnitGcContentFile");

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("getFilenameExtension", "AromaUnitGcContentFile", function(static, ...) {
  "ugc";
}, static=TRUE, protected=TRUE);


setMethodS3("getExtensionPattern", "AromaUnitGcContentFile", function(static, ...) {
  "[.](ugc)$";
}, static=TRUE, protected=TRUE)


setMethodS3("getDefaultColumnNames", "AromaUnitGcContentFile", function(this, ...) {
  "gcContent";
}, protected=TRUE)


setConstructorS3("AromaUgcFile", function(...) {
  this <- extend(AromaUnitGcContentFile(...), "AromaUgcFile");

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


############################################################################
# HISTORY:
# 2009-03-22
# o Created from AromaUflFile.R.
############################################################################
