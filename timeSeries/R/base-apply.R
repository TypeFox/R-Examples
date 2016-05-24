#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  ../../COPYING


################################################################################
# FUNCTION:              DESCRIPTION:
#  apply                  Applies a function to blocks of a 'timeSeries'
################################################################################


setMethod("apply", "timeSeries",
  function(X, MARGIN, FUN, ...)
{
    # A function implemented by Siethelm Wuertz and Yohan Chalabi
    
    # Description:
    #   Apply Functions Over 'Array'timeSeries' Margins
    
    # Arguments:
    #   X  - an array, including a matrix.
    #   MARGIN - a vector giving the subscripts which the function 
    #     will be applied over. E.g., for a matrix 1 indicates rows, 
    #     2 indicates columns, c(1, 2) indicates rows and columns. 
    #     Where X has named dimnames, it can be a character vector 
    #     selecting dimension names.
    #   FUN	- the function to be applied: see ???Details???. In the case 
    #     of functions like +, %*%, etc., the function name must be 
    #     backquoted or quoted.
    #   ...	- optional arguments to FUN.

    # Value:
    #   Returns a vector or array or list of values obtained by 
    #   applying a function to margins of a 'timeSeries'. If the
    #   returned value is a matrix, and if the input argument X and
    #   the returned value have the same number of rows, then the
    #   returned value will be transformed into a 'timeSeries' object.

    # FUNCTION
    
    # Check arguments:
    stopifnot(is.timeSeries(X))
    
    # Extract Title and Documentation:
    Title <- X@title
    Documentation <- X@documentation
    
    # Settings:
    pos <- X@positions
    rec <- X@recordIDs
    FinCenter <- finCenter(X)
    X <- getDataPart(X)
    ans <- callGeneric()
     
    # Manage when univariate timeSeries drops the apply to vector:
    if( is(ans, "vector") && identical(length(ans), NROW(X)) ) 
    {
        ans <- matrix(ans, ncol=1) 
    }
    
    # Result:   
    if (is(ans, "matrix") && identical(NROW(ans), NROW(X))) 
    {
        # Compose timeSeries
        ans <- timeSeries(
          data = ans, charvec = pos,
          one = FinCenter, FinCenter = FinCenter, recordIDs = rec)
        # Preserve Title and Documentation:
        ans@title <- Title
        ans@documentation <- Documentation
    }
    
    # Return Value:
    ans
})


###############################################################################

