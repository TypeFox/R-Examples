setConstructorS3("FileByteVector", function(...) {
  extend(FileVector(..., bytesPerCell=1, storageMode="integer"), "FileByteVector")
})

############################################################################
# HISTORY:
# 2006-01-27
# o Created.
############################################################################ 
