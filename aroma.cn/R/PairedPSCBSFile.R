setConstructorS3("PairedPSCBSFile", function(...) {
  extend(GenericDataFile(...), "PairedPSCBSFile");
})

setMethodS3("loadFit", "PairedPSCBSFile", function(this, ...) {
  pathname <- getPathname(this);
  fit <- loadObject(pathname);

  # Sanity check
  fit <- Arguments$getInstanceOf(fit, "PairedPSCBS");

  fit;
}, protected=TRUE); 


setMethodS3("loadObject", "PairedPSCBSFile", function(this, ...) {
  loadFit(this, ...);
})



#############################################################################
# HISTORY:
# 2021-09-18
# o Added loadFit() == loadObject() for now.
# 2011-01-18
# o Created.
#############################################################################
