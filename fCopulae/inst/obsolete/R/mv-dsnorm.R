
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
# FUNCTION:           DESCRIPTION:
#  dmvsnorm            Multivariate Skew Normal Density Function
#  pmvsnorm            Multivariate Skew Normal Probability Function
#  rmvsnorm            Multivariate Skew Normal Random Deviates
# REQUIREMENTS:       DESCRIPTION:
#  "mvtnorm"           Contributed R - Package
#  "sn" | "mnormt"     Contributed R - Package
################################################################################


################################################################################
# Multivariate Skew Normal Distribution


dmvsnorm = 
function(x, dim = 2, mu = rep(0, dim), Omega = diag(dim), 
alpha = rep(0, dim))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Multivariate Skew Normal Density Function
    
    # Note:
    #   Requires dsn() and dmsn() from R package sn
    
    # FUNCTION:
        
    # Settings:
    xi = mu 
    ans = NA
    
    # Univariate Case:
    if (is.vector(x) & dim == 1) {
        ans = dsn(x, location = xi[1], scale = as.vector(Omega)[1], 
            shape = alpha[1])
    }
    
    # Multivariate Case:
    if (is.matrix(x)) {
        if (dim == ncol(x)) {
            ans = dmsn(x = x, xi = xi, Omega = Omega, alpha = alpha)
        } 
    }
    
    # Check for conflicting Dimensions:
    if (is.na(ans[1])) {
        stop("conflicting x and dim")
    }
        
    # Return Value:
    as.vector(ans)
}


# ------------------------------------------------------------------------------


pmvsnorm = 
function(q, dim = 2, mu = rep(0, dim), Omega = diag(dim), 
alpha = rep(0, dim))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Multivariate Skew Normal Probability Function
    
    # Algorithm:
    
    # Note:
    #   Requires psn() and pmsn() from R package sn
    
    # FUNCTION:
    
    # Settings:
    x = q
    xi = mu
    ans = NA
    
    # Univariate Case:
    if (is.vector(x) & dim == 1) {
        ans = psn(x, location = xi[1], scale = as.vector(Omega)[1], 
            shape = alpha[1])
    }
    
    # Multivariate Case:
    if (is.matrix(x)) {
        if (dim == ncol(x)) {
            ans = NULL
            for (i in 1:nrow(x) ) {
                ans = c(ans, pmsn(x = x[i,], xi = xi, Omega = Omega, 
                    alpha = alpha))
            }
        } 
    }
    
    # Check for conflicting Dimensions:
    if (is.na(ans[1])) {
        stop("conflicting x and dim")
    }
        
    # Return Value:
    as.vector(ans)
}


# ------------------------------------------------------------------------------


rmvsnorm = 
function(n, dim = 2, mu = rep(0, dim), Omega = diag(dim), 
alpha = rep(0, dim))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Multivariate Skew Normal Random Number Generator
    
    # Algorithm:
    
    # Note:
    #   Requires rsn() and rmsn() from R package sn
    
    # FUNCTION:
    
    # Settings:
    ans = NA
    xi = mu
    
    # Univariate Case:
    if (dim == 1) {
        ans = as.matrix(rsn(n, location = xi[1], 
            scale = as.vector(Omega)[1], shape = alpha[1]))
    }
    
    # Multivariate Case:
    if (dim > 1) {
        ans = rmsn(n, xi = xi, Omega = Omega, alpha = alpha)
    }
    
    # Check for conflicting Dimensions:
    if (is.na(ans[1])) {
        stop("dim must be greater 1")
    }
        
    # Return Value:
    rownames(ans) = as.character(1:n)
    colnames(ans) = as.character(1:dim)
    ans
}


################################################################################

