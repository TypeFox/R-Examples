
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

# Copyrights (C)
# for this R-port:
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# METHOD:                 EXTRACTORS:
#  volatility.fGARCH       Returns conditional volatilities for 'fGARCH' objects
################################################################################


volatility.fGARCH <-
    ## better to use S3 style because volatility is defined as a S3 generic
    ## setMethod(f = "volatility", signature(object = "fGARCH"), definition =
    function(object, type = c("sigma", "h"), ...)
{

    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns conditional volatilities for 'fGARCH' objects

    # Arguments:
    #   object - an object of class 'fGarch' as returned by the function
    #       garchFit
    #   type - a character string denoting if the conditional standard
    #       deviations "sigma" or the variances "h" should be returned.
    #   ... - optional argument to be passed, not used.

    # Note:
    #   "volatility" is a generic function. It's default method calculates
    #   (x-mean(x))^2.

    # FUNCTION:

    # Match Arguments:
    type = match.arg(type)

    # Numeric vectors of conditional values:
    if (type == "sigma") {
        volatility = slot(object, "sigma.t")
    } else if (type == "h") {
        volatility = slot(object, "h.t")
    }

    attr(volatility, "type") <- type

    # Return Value:
    volatility

}
##)


################################################################################

