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
#  comment, timeSeries       Get documentation slot of a timeSeries object
#  comment<-,timeSeries      Set documentation slot of a timeSeries object
################################################################################


setMethod("comment", "timeSeries", 
    function(x) 
    {
        # A function implemented by Yohan Chalabi and Diethelm Wuertz

        # Return Value:
        x@documentation
    }
)


# ------------------------------------------------------------------------------

    
setMethod("comment<-", "timeSeries",
    function(x, value)
    {
        x@documentation <- paste(value, collapse = " ")

        # Return Value:
        x
    }
)


################################################################################

