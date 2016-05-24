# Added '...' to some base functions. These will later be
# turned into default functions by setMethodS3().


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Methods in 'base'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# USED TO DO: attach <- appendVarArgs(attach)
attach <- function(...) UseMethod("attach");
setMethodS3("attach", "default", function(...) {
  base::attach(...);
})

# USED TO DO: detach <- appendVarArgs(attach)
detach <- function(...) UseMethod("detach");
setMethodS3("detach", "default", function(...) {
  base::detach(...);
})

# USED TO DO: load <- appendVarArgs(load)
load <- function(...) UseMethod("load");
setMethodS3("load", "default", function(..., envir=parent.frame()) {
  base::load(..., envir=envir);
})

# USED TO DO: save <- appendVarArgs(load)
save <- function(...) UseMethod("save");
setMethodS3("save", "default", function(..., envir=parent.frame()) {
  base::save(..., envir=envir);
})

# USED TO DO: gc <- appendVarArgs(gc)
gc <- function(...) UseMethod("gc");
setMethodS3("gc", "default", function(...) {
  base::gc(...);
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Methods in 'methods'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getClasses <- appendVarArgs(getClasses)
getMethods <- appendVarArgs(getMethods)



############################################################################
# HISTORY:
# 2012-06-20
# o CLEANUP: Dropped non-used adjusted getClass() generic function.
# 2012-02-29
# o Replaced all appendVarArgs() for 'base' functions that do .Internal()
#   calls, because they would then appear as local functions of this
#   package and hence not be accepted by CRAN according to their new
#   policies.  Instead we now create "default" functions that are
#   wrappers to the corresponding functions in the 'base' package.
#   Extra care has to be taken for functions that have arguments whose
#   values are dependent on the call environment/closure.
# o CLEANUP: Dropped unused(?) 'environment <- appendVarArgs(environment)'.
# 2005-02-15
# o Created to please R CMD check.
############################################################################
