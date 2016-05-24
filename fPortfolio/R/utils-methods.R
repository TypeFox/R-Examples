
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


###############################################################################
# FUNCTION:
#  print.solver
#  summary.solver
###############################################################################


print.solver <- 
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Number of Variables:
    nSolution <- length(x$solution)
    
    # Print:
    cat("\nSolver:      ", x$solver)
    cat("\nSolution:    ", 1, ":", x$solution[1])
    for (i in 2:nSolution) 
    cat("\n             ", i, ":", x$solution[i])
    cat("\nObjective:   ", x$objective)
    cat("\nStatus:      ", x$status)
    cat("\nMessage:     ", x$message)
    cat("\n")
}


# -----------------------------------------------------------------------------


.summary.solver <- 
     function(object, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Print:
    print(object[1])
}


###############################################################################

