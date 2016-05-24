
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
#  Ops.timeDate              Group 'Ops' operations on 'timeDate' objects
#  +.timeDate                Performs + operation on 'timeDate' objects
#  -.timeDate                Performs - operation on 'timeDate' objects
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("Ops", c("timeDate", "timeDate"),
    function(e1, e2)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Uses group 'Ops' generic functions for 'timeDate' objects

    # Arguments:
    #   e1 - an object of class 'timeDate'
    #   e2 - an object of class 'timeDate'

    # Value:
    #   Returns the 'Ops' grouped object.

    # FUNCTION:
    ans <- callGeneric(e1@Data, e2@Data)

    if (inherits(ans, "POSIXt"))
        ans <- timeDate(as.character(ans),
                        zone = "GMT", FinCenter = e1@FinCenter)

    # Return Value:
    ans
})


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("+", c("timeDate", "numeric"),
    function(e1, e2)
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz
    
    ans <- callGeneric(e1@Data, e2)
    ans <- timeDate(ans, zone = "GMT", FinCenter = e1@FinCenter)
    
    # Return Value:
    ans
})


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("+", c("numeric", "timeDate"),
    function(e1, e2)
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz
    
    ans <- callGeneric(e1, e2@Data)
    ans <- timeDate(ans, zone = "GMT", FinCenter = e2@FinCenter)

    # Return Value:
    ans
})


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("+", c("timeDate", "timeDate"),
    function(e1, e2) 
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz
    
    stop("binary '+' is not defined for \"timeDate\" objects")
})


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("-", c("timeDate", "numeric"),
    function(e1, e2)
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz
    
    ans <- callGeneric(e1@Data, e2)
    ans <- timeDate(ans, zone = "GMT", FinCenter = e1@FinCenter)
    
    # Return Value:
    ans
})


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("-", c("numeric", "timeDate"),
    function(e1, e2)
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz
    
    stop("Can only subtract from timeDate objects")
})


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("-", c("timeDate", "timeDate"),
    function(e1, e2) 
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz
    
    callGeneric(e1@Data, e2@Data)
})

          
################################################################################
