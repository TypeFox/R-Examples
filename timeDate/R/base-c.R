
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
# MEHODS:                   DESCRIPTION:
#  c.timeDate                Concatenates 'timeDate' objects
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
c.timeDate <-
    function(..., recursive = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Concatenates objects of class 'timeDate'

    # Arguments:
    #   ... - objects to be concatenated.
    #   recursive - a logical. If 'recursive=TRUE', the function
    #       recursively descends through lists combining all their
    #       elements into a vector.

    # Value:
    #   Returns all arguments to be coerced to a 'timeDate' object
    #   which is the type of the returned value.

    # Example:
    #   timeCalendar()[0]; c(timeCalendar()[0], timeCalendar())

    # Details:
    #   This is a generic function which combines its arguments.

    # FUNCTION:

    # List all:
    z = list(...)

    data <- unlist(lapply(z, function(td) c(as.numeric(td, unit = "secs"))))
    timeDate(data, zone = "GMT", FinCenter = z[[1]]@FinCenter)
 }


################################################################################

