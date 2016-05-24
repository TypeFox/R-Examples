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
# FUNCTION:             DESCRIPTION:
#  smoothSupsmu          Smoothes a timeSeries with the supsmu function
#  smoothLowess          Smoothes a timeSeries with the lowess function
#  smoothSpline          Smoothes a timeSeries with the smooth.spline function
# DEPRECATED:           DESCRIPTION:
#  .supsmuSmoother       Smoothes a timeSeries with the supsmu function
#  .lowessSmoother       Smoothes a timeSeries with the lowess function
#  .splineSmoother       Smoothes a timeSeries with the smooth.spline function
################################################################################


# DW: 
# These are older functions which have to be rewritten ...
# The functions are thought to be used to smooth financial 
# price or index series.


# ------------------------------------------------------------------------------


smoothSupsmu <- 
    function(x, bass = 5, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Smoothes a time series with the supsmu function
    
    # Arguments:
    #   x - an univariate timeSeries object, e.g. a price or index series
    #   bass - controls the smoothness of the fitted curve. Values of up 
    #       to 10 indicate increasing smoothness.
    #   ... - further arguments passed to the function supsmu()
    
    # Example:
    #   x <- smoothSupsmu(MSFT[, 4], bass = 0.1); x; plot(x)
    
    # FUNCTION:
    
    # Settings:
    stopifnot(isUnivariate(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
    
    # Handle Missing Values:
    x <- na.omit(x)
    
    # Convert to Vector:
    X <- x
    x <- as.vector(x)
    
    # Smooth:
    ans <- stats::supsmu(x = 1:length(x), y = x, bass = bass, ... )
    data <- cbind(x, ans$y) 
    colnames(data) <- c(colnames(X), "supsmu")
    rownames(data) <- as.character(time(X))
    series(X) <- data
      
    # Preserve Title and Documentation:
    X@title <- Title
    X@documentation <- Documentation
    
    # Return Value:
    X
}


# ------------------------------------------------------------------------------


smoothLowess <- 
    function(x, f = 0.5, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Smoothes a time series with the lowess function
    
    # Arguments:
    #   x - an univariate timeSeries object, e.g. a price or index series
    #   f - the smoother span. This gives the proportion of points in the 
    #       plot which influence the smooth at each value. Larger values 
    #       give more smoothness.
    #   ... - further arguments passed to the function lowess()
    
    # Example:
    #   x = smoothLowess(MSFT[, 4], f = 0.05); x; plot(x)
    
    # FUNCTION:
    
    # Settings:
    stopifnot(isUnivariate(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
    
    # Handle Missing Values:
    x <- na.omit(x)
    
    # Convert to Vector:
    X <- x
    x <- as.vector(x)
    
    # Smooth:
    ans <- stats::lowess(x, f = f, ...)$y
    data <- cbind(x, ans) 
    colnames(data) <- c(colnames(X), "lowess")
    rownames(data) <- as.character(time(X))
    series(X) <- data
      
    # Preserve Title and Documentation:
    X@title <- Title
    X@documentation <- Documentation
    
    # Return Value:
    X
}


# ------------------------------------------------------------------------------


smoothSpline <- 
    function(x, spar = NULL, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Smoothes a time series with the smooth.spline function
    
    # Arguments:
    #   x - an univariate timeSeries object, e.g. a price or index series
    #   f - the smoother span. This gives the proportion of points in the 
    #       plot which influence the smooth at each value. Larger values 
    #       give more smoothness.
    #   ... - further arguments passed to the function smooth.spline()
    
    # Details:
    #   smooth.spline(x, y = NULL, w = NULL, df, spar = NULL, cv = FALSE, 
    #     all.knots = FALSE, nknots = NULL, keep.data = TRUE, df.offset = 0, 
    #     penalty = 1, control.spar = list()) 
   
    # Example:
    #   x = smoothSpline(MSFT[, 4], spar = NULL); x; plot(x)
    #   x = smoothSpline(MSFT[, 4], spar = 0.4); x; plot(x)
    
    # FUNCTION:
    
    # Settings:
    stopifnot(isUnivariate(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
    
    # Handle Missing Values:
    x <- na.omit(x)
    
    # Convert to Vector:
    X <- x
    x <- as.vector(x)
    
    # Smooth:
    ans <- stats::smooth.spline(x, spar = spar, ...)$y
    data <- cbind(x, ans) 
    colnames(data) <- c(colnames(X), "spline")
    rownames(data) <- as.character(time(X))
    series(X) <- data
      
    # Preserve Title and Documentation:
    X@title <- Title
    X@documentation <- Documentation
    
    # Return Value:
    X
}


################################################################################


.supsmuSmoother <-
    function(...)
{
    # FUNCTION:
    
    # Deprecated:
    .Deprecated("smoothSupsmu")
    
    # Return Value:
    smoothSupsmu(...)
}


# ------------------------------------------------------------------------------


.lowessSmoother <-
    function(...)
{
    # FUNCTION:
    
    # Deprecated:
    .Deprecated("smoothLowess")
    
    # Return Value:
    smoothLowess(...)
}


# ------------------------------------------------------------------------------


.splineSmoother <-
    function(...)
{
    # FUNCTION:
    
    # Deprecated:
    .Deprecated("smoothSpline")
    
    # Return Value:
    smoothSpline(...)
}


################################################################################

