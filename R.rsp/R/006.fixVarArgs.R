# Added '...' to some base functions. These will later be
# turned into default functions by setMethodS3().

write <- appendVarArgs(write);

if (exists("restart", mode="function")) {
  restart <- NULL; rm(list="restart"); # To please R CMD check on R (>= 2.15.0)
  restart <- appendVarArgs(restart);
}

# To avoid creating an internal copy with a .Internal() call
stop <- function(...) UseMethod("stop");
setMethodS3("stop", "default", function(...) {
  base::stop(...);
})

############################################################################
# HISTORY:
# 2005-08-01
# o Created to please R CMD check.
############################################################################
