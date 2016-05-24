# This is a hidden function of the l2boost package.

# VAR is a helper function that specifically returns NA if all 
# values of the argument x are NA, otherwise, it returns a var 
# object.
# 
# @param x return variance of x matrix.
#
VAR <- function(x) {if (all(is.na(x))) return(NA) else return(var(x, na.rm=T))}
