###########################################################################/**
# @RdocDefault findAnnotationDataByChipType
#
# @title "Locates an annotation data file by its chip type"
#
# \description{
#  @get "title".
# }
# 
# @synopsis
#
# \arguments{
#   \item{chipType}{A @character string.}
#   \item{pattern}{A filename pattern to search for.}
#   \item{...}{Arguments passed to @see "findAnnotationData".}
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("findAnnotationDataByChipType", "default", function(chipType, pattern=chipType, ...) {
  findAnnotationData(name=chipType, pattern=pattern, set="chipTypes", ...);
}, protected=TRUE)


############################################################################
# HISTORY:
# 2007-09-15
# o Now findAnnotationDataByChipType() utilizes findAnnotationData().
# 2007-02-23
# o BUG FIX: Latest updated of findAnnotationDataByChipType() would not 
#   search recursively.
# 2007-02-21
# o Added Rdoc comments.
# o Added verbose.
# o Added support for aliases.
# o Changed settings$paths$annotationData to settings$annotationData$paths.
# 2007-02-06
# o Created.
############################################################################
