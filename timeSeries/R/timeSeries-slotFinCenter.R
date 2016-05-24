
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
# FUNCTION:                DESCRIPTION:
#  finCenter,timeSeries     Get financial center slot from a 'timeSeries'  
#  finCenter<-,timeSeries   Set financial center slot from a 'timeSeries'  
# FUNCTION:                DESCRIPTION:
#  getFinCenter             Get financial center slot from a 'timeSeries'  
#  setFinCenter<-           Set new financial center slot from a 'timeSeries'  
################################################################################


setMethod("finCenter", "timeSeries",
    function(x)
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz

    # Description:
    #   Get financial center slot from a 'timeSeries' object

    # Arguments:
    #   x - an object of class 'timeSeries'

    # FUNCTION:

    # Extract financial center:
    ans <- x@FinCenter

    # Return Value:
    ans
})


# ------------------------------------------------------------------------------


setMethod("finCenter<-", "timeSeries",
    function(x, value)
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz

    # Description:
    #   Set financial center slot from a 'timeSeries' object

    # Arguments:
    #   x - an object of class 'timeSeries'
    #   value - a character string, setting the name of the financial
    #       center.

    # FUNCTION:
    
    # Check:
    if (x@format == "counts")
        stop(as.character(match.call())[1],
             " is for time series and not for signal series.")

    # Convert to user financial centre:
    positions <- timeDate(charvec = time(x), zone = finCenter(x),
        FinCenter = value)

    # Assign Positions to Time Stamps:
    time(x) <- positions

    # Return Value:
    x
})


################################################################################


getFinCenter <- 
function(x)
{
    # Description:
    #    Get financial center slot from a 'timeSeries' object 
 
    # Arguments:
    #   x - a 'timeSeries' object
    
    # FUNCTION:
    
    # Return Value:
    finCenter(x)
}


# ------------------------------------------------------------------------------


"setFinCenter<-" <- 
function(x, value)
{
    # Description:
    #    Set new financial center slot from a 'timeSeries' object
    
    # FUNCTION:
    
    # Arguments:
    #   x - a 'timeSeries' object
    
    # Assign Financial Center Slot:
    finCenter(x) <- value
    
    # Return Value:
    x
}


################################################################################

