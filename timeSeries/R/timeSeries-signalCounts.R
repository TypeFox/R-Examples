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
#  .signalCounts             Creates charvec for integer indexed time stamps 
################################################################################

   
.signalCounts <-
    function(int)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Descriptions:
    #   Creates the charvec for integer indexed time stamps

    # Arguments:
    #   int - a vector of integers, the counts.

    # FUNCTION:

    # Check that int is an integer
    #   ...

    # Check that all int's are positive ...
    #   ...

    # Format:
    cint <- as.character(int)
    ans <- format(cint, width = max(nchar(cint)), justify = "right")

    # Return Value:
    ans
}


################################################################################

