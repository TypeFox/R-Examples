
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
# FUNCTION:               DESCRIPTION:
#  plot,timeDate           Plots a 'timeDate' object
#  points,timeDate         Adds points to a 'timeDate' plot
#  lines,timeDate          Adds lines to a 'timeDate' plot
#  axis.timeDate           S3 Adds an axis to a 'timeDate' plot
#  pretty.timeDate         S3 Returns a sequence of equally spaced round values 
#  abline,timeDate         Adds an abline to a 'timeDate' plot
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("plot", "timeDate",
    function(x, y, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Plots a 'timeDate' object
    
    # Note:
    #   Doesn't yet support the features of timeDate objects ...

    # FUNCTION:

    # Plot:
    plot(as.POSIXct(x), y, ...)
})


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("points", "timeDate",
    function(x, y, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Note:
    #   Doesn't yet support the features of timeDate objects ...

    # Add Points:
    points(as.POSIXct(x), y, ...)
})


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("lines", "timeDate",
    function(x, y, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds lines to a 'timeDate' plot
    
    # FUNCTION:

    # Note:
    #   Doesn't yet support the features of timeDate objects ...

    # Add Lines:
    lines(as.POSIXct(x), y, ...)
})


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
axis.timeDate <-
    function(side, x, at, format = NULL, labels = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds an axis to a 'timeDate' plot
    
    # Arguments:
    #   side - an integer specifying which side of the plot the axis
    #       is to be drawn on. The axis is placed as follows:
    #       1=below, 2=left, 3=above and 4=right.
    #   x - a 'timeDate' object
    #   at - a 'timeDate' object
    #   format - format string
    #   labels - either a logical value specifying whether annotations
    #       are to be made at the tickmarks, or a vector of character
    #       strings to be placed at the tickpoints.
    #   ... - further arguments to be passed from or to other methods,
    #       typically graphical parameters or arguments of plot.default.
    #       For the plot methods, also format.

    # Notes:
    #   Note that axis.timeDate is not an S3 method !
    #   cannot hence be defined as a S4 method

    # FUNCTION:

    # Format:
    if (is.null(format)) format = whichFormat(x)

    # Add Axis:
    axis.POSIXct(side=side, x=as.POSIXct(x), at=as.POSIXct(at),
        format=format, labels=labels, ...)

    # Return Value:
    invisible()
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
pretty.timeDate <- 
    function(x, n = 5, min.n = n%/%3, shrink.sml = 0.75, 
        high.u.bias = 1.5, u5.bias = 0.5 + 1.5 * high.u.bias, 
        eps.correct = 0, ...) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Returns a sequence of equally spaced round values.
    
    # Details:
    #    Computes a sequence of about n+1 equally spaced ?round? 
    #    values which cover the range of the values in x. 
    #    The values are chosen so that they are 1, 2 or 5 times 
    #    a power of 10.
    
    # Arguments:
    #    x - a 'timeDate' object from which the time is
    #        extracted
    #    n - integer giving the desired number of intervals.  
    #    min.n  - nonnegative integer giving the minimal 
    #        number of intervals. 
    #    shrink.sml - positive numeric by a which a default 
    #        scale is shrunk in the case when range(x) is 
    #        very small.
    #    high.u.bias - non-negative numeric, typically > 1. 
    #        Larger high.u.bias values favor larger units.
    #    u5.bias - non-negative numeric multiplier favoring 
    #        factor 5 over 2.
    #    eps.correct - integer code, one of {0,1,2}. If 
    #       non-0, a correction is made at the boundaries.
    #    ... - further arguments for methods.
       
    # FUNCTION:
    
    x <- as.POSIXct(x)
    ans <- pretty(x, n=n, min.n=min.n, shrink.sml=shrink.sml,
        high.u.bias=high.u.bias, u5.bias=u5.bias, 
        eps.correct=eps.correct, ...)
    
    # Return Value:
    as.timeDate(ans)
}
    

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("abline", signature(v = "timeDate"),
    function(a = NULL, b = NULL, h = NULL, v = NULL, reg = NULL,
        coef = NULL, untf = FALSE, ...)
    {
        # Adds an abline to a 'timeDate' plot
        callGeneric(a = a, b = b, h = h, v = as.POSIXct(v), reg = reg,
            coef = coef, untf = untf, ...)
    }
)


################################################################################
