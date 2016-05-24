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
#  .isOHLC                  Is the series an Open-High-Low-Close series?
#  .isOHLCV                 Is the series an Open-High-Low-Close-Volume series?
################################################################################


# DW:
# I think we need a better method to detect if a series is a OHLC(V) series
# or not. A possible approach would be:
# any High >= Open, Close, Low
# any Low  <= Open, Close, High
# Volume >= 0
# number of columns 4(5)


# -----------------------------------------------------------------------------


.isOHLC <-
    function(object) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Is the series an Open-High-Low-Close series?
    
    # Arguments:
    #   object - an object of class 'timeSeries'
    
    # FUNCTION:
    
    colNames <- c("Open", "High", "Low", "Close")
    if (colnames(object)[1:4] == colNames) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}


# ------------------------------------------------------------------------------ 


.isOHLCV <-
    function(object) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Is the series an Open-High-Low-Close-Volume series?
    
    # Arguments:
    #   object - an object of class 'timeSeries'
    
    # FUNCTION:
    
    colNames <- c("Open", "High", "Low", "Close", "Volume")
    if (colnames(object) == colNames) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}


################################################################################

