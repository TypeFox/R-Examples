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
#  rev,timeSeries            Reverts a 'timeSeries' object in time
################################################################################


.rev.timeSeries <-  function(x) x[NROW(x):1,]


setMethod("rev", "timeSeries", function(x) .rev.timeSeries(x))


# until UseMethod dispatches S4 methods in 'base' functions
rev.timeSeries <- function(x) .rev.timeSeries(x)


################################################################################

