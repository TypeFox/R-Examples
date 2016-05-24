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
# S4 METHOD:                DATABASE ATTACHEMENT:
#  attach,timeSeries         Attaches a 'timeSeries' object to the search path
################################################################################


setMethod("attach", "timeSeries",
function(what, pos = 2, name = deparse(substitute(what)),
    warn.conflicts = TRUE)
{   
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Attaches a 'timeSeries' object
    
    # Details:
    #   The function works in the same way as in the case of a 
    #   data.frame, i.e. the return values are vectors.

    # FUNCTION:

    # Return Value:
    callGeneric(as.data.frame(what), pos, name, warn.conflicts)
})


################################################################################

