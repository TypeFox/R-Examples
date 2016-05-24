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
# FUNCTION:                   DESCRIPTION:
# getDocumentation
# setDocumentation
################################################################################
# FUNCTION:                   MANAGING ATTRIBUTES
#  getAttributes
#  setAttributes<-
# INTERNAL FUNCTION:
#  .appendList
################################################################################


getAttributes <- 
  function (obj) 
{
  # A function implemented by Diethelm Wuertz
    
  # Description:
    
  # FUNCTION:
    
  # Check Argument:
  stopifnot(class(obj) == "timeSeries")
  
  # Extract Attributes:
  ans <- attr(obj@documentation, "Attributes")
    
  # Return Value:
  ans
}


# -----------------------------------------------------------------------------


`setAttributes<-` <- 
  function(obj, value)
{
  # A function implemented by Diethelm Wuertz
    
  # Description:
    
  # Example:
  #   obj <- dummySeries(); getAttributes(obj)
  #   setAttributes(obj) <- list(mat=matrix(1:4, ncol=2)); getAttributes(obj)
  #   getAttributes(obj)$mat[[1]]
    
  # FUNCTION:
    
  # Check Arguments:
  stopifnot(class(obj) == "timeSeries")
  stopifnot(is.list(value))
  stopifnot(length(value) == 1)
  stopifnot(!is.null(value))
  
  # Compose New Attribute:
  name <- names(value)
  names(value) <- NULL
  A <- list(value)
  names(A) <- name
  # print(A)

  # Get Already Existing Attribute
  B <- getAttributes(obj)
  if(is.null(B)) B <- list()
  # print(B)

  # Join Attributes:
  JOINED <- sapply(unique(c(names(A), names(B))), 
    function(x) list(c(A[[x]], B[[x]])))
  # print(JOINED)

  # Assign Attribute:
  attr(obj@documentation, "Attributes") <- JOINED
  
  # Return Value:
  obj 
}


# -----------------------------------------------------------------------------


.appendList <- 
  function (A, B) 
{
  # A function implemented by Diethelm Wuertz
 
  # Description:
  #   Appends list B to list A
    
  # Arguments:
  #   A - first named list element
  #   B - second named list element
  
  # FUNCTION:
    
  # Append list B to list A
  JOINED <- sapply(unique(c(names(A), names(B))), 
    function(x) list(c(A[[x]], B[[x]])))
  
  # Return Value:
  JOINED 
}


###############################################################################

