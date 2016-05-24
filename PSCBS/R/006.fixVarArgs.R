append <- function(...) UseMethod("append");
setMethodS3("append", "default", function(...) {
  base::append(...);
})


############################################################################
# HISTORY:
# 2013-10-14 [HB]
# o ROBUSTNESS: The overriding of append() to become a generic
#   function does now call base::append() in the default, instead
#   of copy the latter.  All this will eventually be removed,
#   when proper support for c, [, [[ etc. has been added everywhere.
# 2010-10-02
# o Created to please R CMD check.
############################################################################
