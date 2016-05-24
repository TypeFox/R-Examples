setMethodS3("getChecksum", "Object", function(object, ...) {
  object <- clearCache(object);
  NextMethod("getChecksum", object=object);
}, export=FALSE)


############################################################################
# HISTORY:
# 2014-02-03
# o ROBUSTNESS, now passing the first/dispatch argument (named 'object')
#   to argument 'object' of NextMethod() as a named argument, which just
#   happens to be of the same names.
# 2011-04-02
# o Added getChecksum() for the Object class.
# o Created.
############################################################################
