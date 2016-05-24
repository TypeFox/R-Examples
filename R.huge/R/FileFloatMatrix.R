setConstructorS3("FileFloatMatrix", function(...) {
  extend(FileMatrix(..., bytesPerCell=4, storageMode="double"), "FileFloatMatrix")
})


############################################################################
# HISTORY:
# 2006-02-18
# o Created.
############################################################################ 
