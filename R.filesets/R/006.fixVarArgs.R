append <- function(...) UseMethod("append");
setMethodS3("append", "default", function(...) {
  base::append(...);
})

readLines <- function(...) UseMethod("readLines");
setMethodS3("readLines", "default", function(...) {
  base::readLines(...);
})


############################################################################
# HISTORY:
# 2013-10-07 [HB]
# o ROBUSTNESS: The overriding of append() to become a generic
#   function does now call base::append() in the default, instead
#   of copy the latter.  All this will eventually be removed,
#   when proper support for c, [, [[ etc. has been added everywhere.
# 2012-03-06 [HB]
# o CRAN POLICY: Removed all internal copies of 'base' functions that
#   have .Internal() calls.
# 2007-02-27 [HB]
# o BUG FIX: Removed explicit reference to 'base' etc again. The reason is
#   that if a previous package already modified, say, write(), to become a
#   generic function, that was overwritten again when this package was
#   loaded.
# 2007-02-23 [KS]
# o Make explicit reference to 'base' - this is safer, in case of colMeans()
#   defined by user or other packages.
# 2006-03-24
# o Created to please R CMD check.
############################################################################
