
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
#  skewness                  Returns a number which is the skewness of the data
#  skewness.default          Default method
#  skewness.data.frame       Method for objects of class data.frame
#  skewness.POSIXct          Method for objects of class POSIXct
#  skewness.POSIXlt          Method for objects of class POSIXlt
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
skewness <-
    function (x, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Return Value:
    UseMethod("skewness")
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
skewness.default <-
    function (x, na.rm = FALSE, method = c("moment", "fisher"), ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns the value of the skewness of a distribution function.

    # Details:
    #   Missing values can be handled.

    # FUNCTION:

    # Method:
    method = match.arg(method)

    # Warnings:
    if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
        warning("argument is not numeric or logical: returning NA")
        return(as.numeric(NA))}

    stopifnot(NCOL(x) == 1)

    # Remove NAs:
    if (na.rm) x = x[!is.na(x)]

    # Skewness:
    n = length(x)
    if (is.integer(x)) x = as.numeric(x)

    # Selected Method:
    if (method == "moment") {
        skewness = sum((x-mean(x))^3/sqrt(as.numeric(var(x)))^3)/length(x)
    }
    if (method == "fisher") {
        if (n < 3)
            skewness = NA
        else
            skewness = ((sqrt(n*(n-1))/(n-2))*(sum(x^3)/n))/((sum(x^2)/n)^(3/2))
    }

    # Add Control Attribute:
    attr(skewness, "method") <- method

    # Return Value:
    skewness
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
skewness.data.frame <-
    function (x, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Return Value:
    sapply(x, skewness, ...)
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
skewness.POSIXct <-
    function (x, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Return Value:
    structure(skewness(unclass(x), ...), oldClass(x))
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
skewness.POSIXlt <-
    function (x, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Return Value:
    as.POSIXlt(skewness(as.POSIXct(x), ...))
}


################################################################################

