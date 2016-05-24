######################################################################
#
# assignment.R
#
# Written by Carter T. Butts <buttsc@uci.edu>; portions contributed by
# David Hunter <dhunter@stat.psu.edu> and Mark S. Handcock
# <handcock@u.washington.edu>.
#
# Last Modified 05/04/14
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/network package
#
# This file contains various routines for the assignment of network objects
# into calling environments. These are internal functions and not to be used
# by the package users.
#
# Contents:
#
#   .findNameInSubsetExpr
#   .validLHS
#
######################################################################


# Recursively traverse the parse tree of the expression x, ensuring that it is
# a valid subset expresssion, and return the name associated with the expression.
#
.findNameInSubsetExpr <- function(x){
  if (class(x)=='call'){
    # Ensure call is a subset function, one of $, [, or [[
    if(!(deparse(x[[1]]) %in% c('$','[','[['))) return(NA)

    # Make sure arguments are clean
    xns <- lapply(x[2:length(x)],.findNameInSubsetExpr)
    if (any(is.na(xns))) return(NA)
      
    # Possible name found
    return(xns[[1]])
  }
  else if (class(x)=='name')
    return(deparse(x))
 
  NULL 
}

# Return TRUE if x is a valid left-hand-side object that can take a value

.validLHS <- function(x,ev){
  xn <- .findNameInSubsetExpr(x)
  # There are valid expressions for which we don't want to assign into the caller's env.
  # For instance, when a user executes z<-add.edges(x+y), then the user obviously
  # doesn't want x+y to be assigned. Rather he's using them as temporaries to obtain
  # z. OTOH we don't want someone doing something obtuse like add.edges(x[sample(...)])
  # In the first case, it's not wrong to end up here, but in the second case we would
  # like to warn the user. But we're not going to at this point.
  #warning('Cannot make assignment into ',deparse(x))
  if (!is.null(xn) && !is.na(xn) && exists(xn,envir=ev))
    return(TRUE)
  else
    return(FALSE)
}
