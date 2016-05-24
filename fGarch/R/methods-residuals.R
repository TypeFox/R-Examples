
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
#  residuals.fGARCH        S4 residuals method for an object of class 'fGARCH'
################################################################################


setMethod(f = "residuals", signature(object = "fGARCH"), definition =
    function(object, standardize = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   S4 Residuals method for an object of class fGARCH

    # Arguments:
    #   object - an object of class fGarch as returned by the
    #       function garchFit
    #   ... - optional argument to be passed, this may be standardize=FALSE
    #       to return the -non-standardized values of the residuals.
    #       By default standardized residuals will be returned.

    # FUNCTION:


    # Residuals:
    if (standardize) {
        ans = object@residuals/object@sigma.t
    } else {
        ans = object@residuals
    }

    # Return Value:
    ans

})


# ------------------------------------------------------------------------------
