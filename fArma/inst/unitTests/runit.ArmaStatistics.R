
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port: 
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:               TRUE ARMA STATISTICS:
#  armaTrueacf             Returns True ARMA autocorrelation function
#  armaRoots               Returns Roots of the ARMA characteristic polynomial
################################################################################


test.armaTrueacf = 
function()
{ 
    # armaTrueacf: Returns True ARMA autocorrelation function

    # armaTrueacf(model, lag.max = 20, type = "correlation", doplot = TRUE)
    model = list(ar = c(0.5, -0.5))
    armaTrueacf(model, lag.max = 10)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.armaRoots = 
function()
{ 
    # armaRoots:   Returns Roots of the ARMA characteristic polynomial
    
    # armaRoots(coefficients, n.plot = 400, digits = 4, ...)
    coefficients = c(-0.5, 0.9, -0.1, -0.5)
    ans = armaRoots(coefficients)
    target = round(sum(ans), 2)
    checkSum = 4.58
    checkEqualsNumeric(target, checkSum)
    
    # Return Value:
    return()    
}


################################################################################

