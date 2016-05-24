setMethodS3("getDeviceResolution", "default", function(scale=1, ...) {
  res <- scale * par("cra") / par("cin");
  res;
})


############################################################################
# HISTORY:
# 2007-01-25
# o Created.
############################################################################
