
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
#  start.timeDate            Extracts the first entry of a 'timeDate' object
#  end.timeDate              Extracts the last entry of a 'timeDate' object
#  min.timeDate              Extracts the smallest entry of a 'timeDate' object
#  max.timeDate              Extracts the largest entry of a 'timeDate' object
#  range.timeDate            Extracts range of a 'timeDate' object
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
start.timeDate <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Extracts the first object of a 'timeDate' vector

    # Arguments:
    #   x - a 'timeDate' object

    # Value:
    #   Returns from 'x' the earliest entry as an object of class
    #   'timeDate'.

    # FUNCTION:

    # Return Value:
    timeDate(min(as.POSIXct(x)), zone = "GMT", FinCenter = x@FinCenter)
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
end.timeDate <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Extracts the last object of a 'timeDate' vector

    # Arguments:
    #   x - a 'timeDate' object

    # Value:
    #   Returns an object of class 'timeDate'.

    # FUNCTION:

    # Return Value:
    timeDate(max(as.POSIXct(x)), zone = "GMT", FinCenter = x@FinCenter)
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
min.timeDate <- 
    function(..., na.rm = FALSE) 
{
    # DW: is this ok ??
    
    start.timeDate(...)
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
max.timeDate <- 
    function(..., na.rm = FALSE) 
{
    # DW: is this ok ??
    
    end.timeDate(...)
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
range.timeDate <- 
    function(..., na.rm = FALSE) 
{
    c(start(..., na.rm = na.rm), end(..., na.rm = na.rm)) 
}


################################################################################

