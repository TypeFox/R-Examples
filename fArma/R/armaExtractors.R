
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
# FUNCTION:               SIMULATION AND FITTING:
#  coef.fARMA              S3: Returns coefficidents from a fitted ARMA object
#  coefficients.fARMA      S3: Synonyme for coef.fARMA
#  fitted.fARMA            S3: Returns fitted values from a fitted ARMA object
#  residuals.fARMA         S3: Returns residuals from a fitted ARMA object
################################################################################


coef.fARMA =
function(object, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns coefficients from a fitted ARMA object
    
    # Note:
    #   Alternatively you can use coefficient().

    # FUNCTION:
    
    # Coefficients:
    ans = object@fit$coef
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


fitted.fARMA = 
function(object, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns fitted values from a fitted ARMA object

    # FUNCTION:
        
    # Fitted Values:
    ans = object@fitted$fitted
    classAns = class(object@data$x)
    if (classAns == "ts") ans = as.ts(ans)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


residuals.fARMA = 
function(object, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns Residuals from a Fitted ARMA Object
       
    # FUNCTION:
    
    # Check:
       
    # Residual Values:
    ans = object@residuals$residuals
    classAns = class(object@data$x)
    if (classAns == "ts") ans = as.ts(ans)
    
    # Return Value:
    ans
}


################################################################################

