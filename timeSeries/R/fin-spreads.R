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
#  spreads                   Computes spreads from a 'timeSeries' object
#  midquotes                 Computes mid quotes from a 'timeSeries' object
################################################################################


# DW:
# Setting bid and ask for column names is maybe the best choice. Examples
# are the TED spread or the Libo OIS spread. The spread between High and Low
# is the range.


# ------------------------------------------------------------------------------


spreads <-
    function(x, which = c("Bid", "Ask"), tickSize = NULL)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes spreads from a 'timeSeries' object
     
    # FUNCTION:
      
    # Check arguments:
    stopifnot(is.timeSeries(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
    
    # Compute Spread:
    spread <- x[, which[2]] - x[, which[1]]
    if (!is.null(tickSize)) series(spread) <- round(series(spread)/tickSize)

    # Preserve Title and Documentation:
    spread@title <- Title
    spread@documentation <- Documentation
      
    # Return Value:
    spread
}


# ------------------------------------------------------------------------------


midquotes =
function(x, which = c("Bid", "Ask"))
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes mid quotes from a 'timeSeries' object
     
    # FUNCTION:
    
    # Check arguments:
    stopifnot(is.timeSeries(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
  
    # Compute Mid Quotes:
    midquotes = 0.5 * ( x[, which[1]] + x[, which[2]] )
   
    # Preserve Title and Documentation:
    midquotes@title <- Title
    midquotes@documentation <- Documentation  
  
    # Return Value:
    midquotes
}


################################################################################

