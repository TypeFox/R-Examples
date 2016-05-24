
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
# FUNCTION:                 DESCRIPTION:
#  VaR                       Computes Value-at-Risk
#  CVaR                      Computes Conditional Value-at-Risk
################################################################################


VaR =
function(x, alpha = 0.05, type = "sample", tail = c("lower", "upper"))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Value-at-Risk
    
    # Arguments:
    #   x -  an uni- or multivariate timeSeries object
    #   alpha - a numeric value, the confidence interval
    #   type - a character string, the type to calculate the value-at-risk
    #   tail - a character string denoting which tail will be
    #       considered, either \code{"lower"} or \code{"upper"}.
    #       If \code{tail="lower"}, then alpha will be converted to
    #       \code{alpha=1-alpha}.
    
    # FUNCTION:
    
    # Settings:
    x = as.matrix(x)
    tail = match.arg(tail)
    
    # Value-at-Risk:
    if (type == "sample") {
        if (tail == "upper") alpha = 1-alpha
        # Important: use type=1 !
        VaR = quantile(x, probs = alpha, type = 1)
    } else if (type == "gpd") {
        VaR = "Not yet Implemented"
    } else if (type == "obre") {
        VaR = "Not yet Implemented"
    }
    
    # Return Value:
    VaR
}
   

# ------------------------------------------------------------------------------


CVaR = 
function(x, alpha = 0.05, type = "sample", tail = c("lower", "upper"))  
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Conditional Value-at-Risk
    
    # Arguments:
    #   x -  an uni- or multivariate timeSeries object
    #   alpha - a numeric value, the confidence interval
    #   type - a character string, the type to calculate the value-at-risk
    #   tail - a character string denoting which tail will be considered, 
    #       either "lower" or upper", if tail="lower", then alpha will be 
    #       converted to alpha=1-alpha.
    
    # FUNCTION:
        
    # Settings:
    x = as.matrix(x)
    tail = match.arg(tail)
    
    # Sample VaR:
    VaR = VaR(x, alpha, type, tail)
    
    # Sample CVaR:
    if (tail == "upper") alpha = 1-alpha
    if (type == "sample") { 
         CVaR = NULL
         for (i in 1:ncol(x)) {
            X = as.vector(x[, i])
            CVaR = c(CVaR, 
                VaR[i] - 0.5 * mean(((VaR[i]-X) + abs(VaR[i]-X))) / alpha ) 
        }
    }
    
    # Return Value:
    CVaR
}
    

################################################################################

