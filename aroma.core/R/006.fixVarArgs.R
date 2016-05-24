# Added '...' to some base functions. These will later be
# turned into default functions by setMethodS3().

# USED TO DO: colSums <- appendVarArgs(colSums);
colSums <- function(...) UseMethod("colSums");
setMethodS3("colSums", "default", function(...) {
  base::colSums(...);
})

# USED TO DO: colMeans <- appendVarArgs(colMeans);
colMeans <- function(...) UseMethod("colMeans");
setMethodS3("colMeans", "default", function(...) {
  base::colMeans(...);
})

write <- appendVarArgs(write);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Add generic process().
# NB: process() is defined in R.rsp (>= 0.9.1)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setGenericS3("process", overwrite=TRUE);


############################################################################
# HISTORY:
# 2013-10-07 [HB]
# o CLEANUP: No longer copying matrixStats::colMedians() as a workaround.
# 2013-08-03 [HB]
# o Now aroma.core use the same type of generic colMedians() function
#   as matrixStats.
# 2012-07-20 [HB]
# o Added '...' to write(), because we no longer attach R.rsp.
# 2012-03-01 [HB]
# o Replaced all appendVarArgs() for 'base' functions that do .Internal()
#   calls, because they would then appear as local functions of this
#   package and hence not be accepted by CRAN according to their new
#   policies.  Instead we now create "default" functions that are
#   wrappers to the corresponding functions in the 'base' package.
#   Extra care has to be taken for functions that have arguments whose
#   values are dependent on the call environment/closure.
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
