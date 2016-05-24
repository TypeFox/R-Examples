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
#  scale,timeSeries          Centers and/or scales a 'timeSeries' object
################################################################################


.scale.timeSeries <- 
  function(x, center = TRUE, scale = TRUE)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Centers and/or scales a 'timeSeries' object.

    # Arguments:

    # FUNCTION:

    # Scale:
    setDataPart(x, scale(x = getDataPart(x), center = center, scale = scale))
}


setMethod("scale", "timeSeries",
          function(x, center = TRUE, scale = TRUE)
          .scale.timeSeries(x, center = center, scale = scale))


# until UseMethod dispatches S4 methods in 'base' functions
scale.timeSeries <- function (x, center = TRUE, scale = TRUE)
    .scale.timeSeries(x, center = center, scale = scale)


################################################################################

