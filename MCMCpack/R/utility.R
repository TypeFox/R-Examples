##########################################################################
## Utility Functions
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


# takes a symmetric matrix x and returns lower diagonal
# note: does not check for symmetry
#
# ADM 4/18/2003 

"vech" <-
  function (x) {
    x <- as.matrix(x)
    if (dim(x)[1] != dim(x)[2]) {
      stop("Non-square matrix passed to vech().\n")
    }
    output <- x[lower.tri(x, diag = TRUE)]
    dim(output) <- NULL
    return(output)
  }

# takes vector x and returns an nrow times nrow symmetric matrix
# this will recycle the elements of x as needed to fill the matrix
#
# ADM 4/18/2003
# ADM 11/13/2003 [bug fix]
# ADM 1/25/2006 [patch to automatically compute matrix size]

"xpnd" <-
  function (x, nrow = NULL) {
    dim(x) <- NULL
    if(is.null(nrow)) nrow <- (-1 + sqrt(1 + 8 * length(x))) / 2
    output <- matrix(0, nrow, nrow)
    output[lower.tri(output, diag = TRUE)] <- x
    hold <- output
    hold[upper.tri(hold, diag=TRUE)] <- 0
    output <- output + t(hold)    
    return(output)
  }
