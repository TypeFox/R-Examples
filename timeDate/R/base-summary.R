
# This R package is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This R package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this R package; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                 DESCRIPTION:
#  summary.timeDate          Summarizes details of a 'timeDate' object
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
summary.timeDate <-
    function(object, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Summarizes details of a 'timeDate' object

    # Arguments:
    #   x - a 'timeDate' object to be summarized.

    # Effect:
    #   Produce a summary report of the details of a 'timeDate'
    #   object.

    # Print:
    x = object
    cat(  "Object:       ", as.character(match.call())[2])
    cat("\nStart Record: ", as.character(start(x)))
    cat("\nEnd Record:   ", as.character(end(x)))
    cat("\nObservations: ", length(as.character(x)))
    cat("\nFormat:       ", x@format)
    cat("\nFinCenter:    ", x@FinCenter)
    cat("\n")

    # Return Value:
    invisible(object)
}


################################################################################

