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
#  isUnivariate              Tests if a 'timeSeries' object is univariate
#  isMultivariate            Tests if a 'timeSeries' object is multivariate
################################################################################


isUnivariate <-
  function(x)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Tests if a time series or rectangular object is univariate

    # FUNCTION:

    # Return Value:
    if (NCOL(x) == 1) TRUE else FALSE
}


# ------------------------------------------------------------------------------


isMultivariate <-
function(x)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Tests if a time series or rectangular object is multivariate

    # FUNCTION:

    # Return Value:
    if (NCOL(x) > 1) TRUE else FALSE
}


################################################################################

