# @author "HB"
setConstructorS3("AffymetrixCsvFile", function(..., sep=",", .verify=TRUE) {
  this <- extend(AffymetrixTabularFile(..., .verify=FALSE), "AffymetrixCsvFile");

  if (.verify) {
    verify(this, ...);
  }

  this;
})


setMethodS3("getDefaultExtension", "AffymetrixCsvFile", function(static, ...) {
  "csv";
}, static=TRUE, protected=TRUE)

setMethodS3("getExtensionPattern", "AffymetrixCsvFile", function(static, ...) {
  ext <- getDefaultExtension(static, ...);
  pattern <- sprintf("[.](%s|%s)$", tolower(ext), toupper(ext));
  pattern;
}, static=TRUE, protected=TRUE)


setMethodS3("findByChipType", "AffymetrixCsvFile", function(static, chipType, pattern=sprintf("^%s.*%s", chipType, getExtensionPattern(static)), ...) {
  NextMethod("findByChipType", chipType=chipType, pattern=pattern);
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2011-11-19
# o Added getDefaultExtension() to AffymetrixCsvFile.
# 2007-09-10
# o Created from AffymetrixCsvGenomeInformation.R.
############################################################################
