
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
#  outlier,timeSeries        Removes outliers from a 'timeSeries' object
################################################################################


# DW:
# We should call this function no longer outlier, much better woud be
# splits() since the function tries to detect splits by large outliers.
# For outlier detection we should use better methods than just the sd().


# ------------------------------------------------------------------------------


splits <- 
  function(x, sd = 3, complement = TRUE, ...)
{
    # Return Value:
    outlier(x=x, sd=sd, complement=complement, ...)
}


# ------------------------------------------------------------------------------


setMethod("outlier", "ANY",
    function(x, sd = 3, complement = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns outlier splits

    # Arguments:
    #   x - a numeric vector
    #   sd - a numeric value of standard deviations, e.g. 5
    #       means that values larger or smaller tahn five
    #       times the standard deviation of the series will
    #       be detected.
    #   complement - a logical flag, should the outlier series
    #       or its complements be returned.

    # Note:
    #   This function is thought to find splits in financial
    #   price or index series If a price or index is splitted we
    #   observe in the returns a big jump of several standard
    #   deviations which is identified usual as an outlier.

    # FUNCTION:

    # Check arguments:
    stopifnot(is.timeSeries(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
      
    # Find Outliers:
    SD <- sd * sd(x)
    if (complement) {
        ans  <- x[x <= SD]
    } else {
        ans <- x[x > SD]
        names(ans) <- as.character(which(x > SD))
    }
      
    # Preserve Title and Documentation:
    ans@title <- Title
    ans@documentation <- Documentation

    # Return Value:
    ans
})


# ------------------------------------------------------------------------------


setMethod("outlier", "timeSeries",
    function(x, sd = 3, complement = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns outliers in a timeSeries object or the complement

    # Arguments:
    #   x - an object of class 'timeSeries'.
    #   sd - a numeric value of standard deviations, e.g. 5
    #       means that values larger or smaller tahn ten
    #       times the standard deviation of the series will
    #       be removed.
    #   complement - a logical flag, should the outler series
    #       or its complement be returned.

    # FUNCTION:

    # Check if univariate Series:
    if (!isUnivariate(x))
        stop("Supports only univariate timeSeries Objects")

    # Find Outliers:
    SD = sd * sd(x)
    if (complement) {
        x  = x[abs(x) <= SD,]
    } else {
        x = x[abs(x) > SD,]
    }

    # Return Value:
    x
})


################################################################################

