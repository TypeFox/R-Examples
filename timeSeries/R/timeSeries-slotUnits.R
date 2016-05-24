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
# FUNCTION:                 DESCRIPTION:
#  getUnits                  Get units slot from a 'timeSeries'  
#  setUnits<-                Set new units slot to a 'timeSeries'
################################################################################


getUnits <-
  function(x)
  {
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION: 
    
    # Return Value:
    UseMethod("getUnits")
  }

getUnits.default <- 
  function(x)
  {
    # Description:
    #   Get units slot from a 'timeSeries' object.
    
    # Arguments:
    #   x - a 'timeSeries' object
    
    # FUNCTION:
    
    # Return Value:
    colnames(x)
  }


# ------------------------------------------------------------------------------


"setUnits<-" <-
  function(x, value)
  {
    # Description:
    #   Set units slot to a 'timeSeries' object.
    
    # Arguments:
    #   x - a 'timeSeries' object
    
    # FUNCTION:
    
    # Assign Time Slot:
    colnames(x) <- value
    
    # Return Value:
    x    
  }



################################################################################

