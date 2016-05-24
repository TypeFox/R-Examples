
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
#  sample.timeDate           Resamples a 'timeDate' object
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("sample", "timeDate",
    function(x, size, replace = FALSE, prob = NULL)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # FUNCTION:

    # Sample:
    ct <- sample(as.POSIXct(x), size, replace, prob)
    ans <- timeDate(ct, zone = "GMT", FinCenter = x@FinCenter)

    # Return Value:
    ans
})


################################################################################

