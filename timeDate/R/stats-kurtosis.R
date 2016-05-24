
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
#  kurtosis                  Returns a number which is the kurtosis of the data
#  kurtosis.default          Default method
#  kurtosis.data.frame       Method for objects of class data.frame
#  kurtosis.POSIXct          Method for objects of class POSIXct
#  kurtosis.POSIXlt          Method for objects of class POSIXlt
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
kurtosis <-
    function (x, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Return Value:
    UseMethod("kurtosis")
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
kurtosis.default <-
    function (x, na.rm = FALSE, method = c("excess", "moment", "fisher"), ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns the value of the kurtosis of a distribution function.

    # Details:
    #   Missing values can be handled.

    # FUNCTION:

    # Method:
    method = match.arg(method)

    stopifnot(NCOL(x) == 1)

    # Warnings:
    if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
        warning("argument is not numeric or logical: returning NA")
        return(as.numeric(NA))}

    # Remove NAs:
    if (na.rm) x = x[!is.na(x)]

    # Kurtosis:
    n = length(x)
    if (is.integer(x)) x = as.numeric(x)
    if (method == "excess") {
        kurtosis = sum((x-mean(x))^4/as.numeric(var(x))^2)/length(x) - 3
    }
    if (method == "moment") {
        kurtosis = sum((x-mean(x))^4/as.numeric(var(x))^2)/length(x)
    }
    if (method == "fisher") {
        kurtosis = ((n+1)*(n-1)*((sum(x^4)/n)/(sum(x^2)/n)^2 -
            (3*(n-1))/(n+1)))/((n-2)*(n-3))
    }

    # Add Control Attribute:
    attr(kurtosis, "method") <- method

    # Return Value:
    kurtosis
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
kurtosis.data.frame <-
    function (x, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Return Value:
    sapply(x, kurtosis, ...)
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
kurtosis.POSIXct <-
    function (x, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Return Value:
    structure(kurtosis(unclass(x), ...), class = oldClass(x))
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
kurtosis.POSIXlt <-
    function (x, ...)
{
    # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Return Value:
    as.POSIXlt(kurtosis(as.POSIXct(x), ...))
}


################################################################################

