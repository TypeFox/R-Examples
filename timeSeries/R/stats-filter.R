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
#  filter,timeSeries         Applies linear filtering to a 'timeSeries' object 
################################################################################


setMethod("filter", "timeSeries",
    function(x, filter, method = c("convolution", "recursive"),
        sides = 2, circular = FALSE, init = NULL)
{
    # Description:
    #   Applies linear filtering to a 'timeSeries' object 
    
    # Arguments:
    #   x - a univariate or multivariate time series.
    #   filter - a vector of filter coefficients in reverse time order (as  
    #       for AR or MA coefficients).
    #   method - Either "convolution" or "recursive" (and can be 
    #       abbreviated). If "convolution" a moving average is used: 
    #       if "recursive" an autoregression is used.
    #   sides - for convolution filters only. If sides=1 the filter 
    #       coefficients are for past values only; if sides=2 they are 
    #       centred around lag 0. In this case the length of the filter 
    #       should be odd, but if it is even, more of the filter is 
    #       forward in time than backward.
    #   circular - for convolution filters only. If TRUE, wrap the filter 
    #       around the ends of the series, otherwise assume external 
    #       values are missing (NA).
    #   init - for recursive filters only. Specifies the initial values 
    #       of the time series just prior to the start value, in reverse 
    #       time order. The default is a set of zeros.
    
    # Value:
    #    Returns a 'timeSeries' object.
    
    # FUNCTION:
    
    # Check Arguments:
    stopifnot(is.timeSeries(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
      
    # Filter:
    ans <- filter(getDataPart(x), filter = filter, method = method,
        sides = sides, circular = circular, init = init)
    # Note: do not use as.matrix because ts objects might
    #       not be coerced properly
    ans <- as(ans, "matrix")
    
    # Add Column Names:
    colnames(ans) <- colnames(x)
    ans <- setDataPart(x, ans)
      
    # Preserve Title and Documentation:
    ans@title <- Title
    ans@documentation <- Documentation
      
    # Return Value:
    ans
    
})


################################################################################

