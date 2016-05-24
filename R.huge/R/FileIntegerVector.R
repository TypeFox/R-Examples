setConstructorS3("FileIntegerVector", function(...) {
  extend(FileVector(..., bytesPerCell=4, storageMode="integer"), "FileIntegerVector")
})

############################################################################
# HISTORY:
# 2006-01-27
# o Created.
############################################################################ 
