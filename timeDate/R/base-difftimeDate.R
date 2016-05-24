
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
#  difftimeDate              Returns a difference of two 'timeDate' objects
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
difftimeDate <-
    function(time1, time2,
    units = c("auto", "secs", "mins", "hours", "days", "weeks"))
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Takes a difference of two 'timeDate' objects

    # Arguments:
    #   time1, time2 - 'timeDate' objects.
    #   units - a character string. The units in which the results
    #       are desired, one of the following: "auto", "secs",
    #       "mins", "hours", "days", "weeks"

    # Value:
    #   'difftimeDate' takes a difference of two 'timeDate'
    #   objects and returns an object of class 'difftime' with
    #   an attribute indicating the units.

    # FUNCTION:

    # Convert to GMT:
    time1GMT = timeDate(time1, zone = time1@FinCenter,
        FinCenter = "GMT")
    time2GMT = timeDate(time2, zone = time2@FinCenter,
        FinCenter = "GMT")

    # Return Value:
    difftime(time1GMT@Data, time2GMT@Data, tz = "GMT", units = units[1])
}



################################################################################

