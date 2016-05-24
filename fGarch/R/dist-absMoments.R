
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
# FUNCTION:              MOMENTS:
#  absMoments             Compute absolute moments of a symmetric distribution
################################################################################


absMoments <- 
function(n, density = c("dnorm", "dged", "dstd"), ...)
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Compute the absolute moments of a standardized
    #   symmetric distribution function.
    
    # Arguments:
    #   n - a vector of integers i, to compute M_i
    #   density - a character denoting the density
    #       "norm", "ged", "std" or any other
    #   ... - parameters passed to the standardized
    #       symmetric density function
    
    # Value:
    #   Returns a numeric vector of moments M_i.
    #   Stores globally errors in the variable absMoment.error
    #     if the moments were comuted numerically.
    
    # FUNCTION:
              
    # norm - Normal Distribution:
    if (density == "dnorm" | density == "norm") {
        return (sqrt(2)^n * gamma((n+1)/2) / sqrt(pi)) }

    # ged - Generalized Error Distribution:
    if (density == "dged" | density == "ged") {
        parm = function(n, nu) {
            lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
            return ((2^(1/nu)*lambda)^n * gamma((n+1)/nu) / gamma(1/nu)) }
        return(parm(n, ...)) 
    }

    # std - Student-t Distribution:
    # Note: nu > 2*n
    if (density == "dstd" | density == "std") {
        parm = function(n, nu) {
            return (beta(1/2+2*n, nu/2-2*n)/beta(1/2, nu/2) * sqrt(nu-2)) }
        return(parm(n, ...)) 
    }

    # Any other standardized symmetric Distribution ...
    fun = match.fun(density)
    moments = function(x, n, ...) { 2 * x^n * fun(x, ...) }
    M = .absMoments.error <- NULL
    for (i in n) {
        I = integrate(moments, 0, Inf, n = i, ...)
        M = c(M, I$value)
        .absMoments.error <- c(.absMoments.error, I$abs.error) 
    }
    attr(M, "control") <- .absMoments.error 
    return(M)

    # Return Value:
    invisible()
}


################################################################################

