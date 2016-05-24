setConstructorS3("CbsSegmentationDataFile", function(...) {
  extend(SegmentationDataFile(...), "CbsSegmentationDataFile");
})

setMethodS3("loadFit", "CbsSegmentationDataFile", function(this, ...) {
  pathname <- getPathname(this);
  loadObject(pathname);
}, protected=TRUE);



############################################################################
# HISTORY:
# 2010-08-05
# o Created.
############################################################################
