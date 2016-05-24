setConstructorS3("FileIntegerMatrix", function(...) {
  extend(FileMatrix(..., bytesPerCell=4, storageMode="integer"), "FileIntegerMatrix")
})

############################################################################
# HISTORY:
# 2006-01-22
# o Created.
############################################################################ 
