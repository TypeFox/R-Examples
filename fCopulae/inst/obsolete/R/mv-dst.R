
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
#  dmvst               Multivariate Skew Sudent-t Density Function
#  pmvst               Multivariate Skew Sudent-t Probability Function
#  rmvst               Multivariate Skew Sudent-t Random Deviates
# REQUIREMENTS:       DESCRIPTION:
#  "mvtnorm"           Contributed R - Package
#  "sn" | "mnormt"     Contributed R - Package
################################################################################


################################################################################


dmvst = 
function(x, dim = 2, mu = rep(0, dim), Omega = diag(dim), 
alpha = rep(0, dim), df = 4)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Multivariate Skew Sudent-t Density Function
    
    # Arguments:
    
    # FUNCTION:
    
    # Settings:
    xi = mu
    ans = NA
    
    # Univariate Case:
    if (is.vector(x) & dim == 1) {
        ans = dst(x, location = xi[1], scale = as.vector(Omega)[1], 
            shape = alpha[1], df = Inf)
    }
    
    # Multivariate Case:
    if (is.matrix(x)) {
        if (dim == ncol(x)) {
            ans = dmst(x = x, xi = xi, Omega = Omega, alpha = alpha, df = df)
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


pmvst = 
function(q, dim = 2, mu = rep(0, dim), Omega = diag(dim), 
alpha = rep(0, dim), df = 4)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Multivariate Skew Sudent-t Probability Function
    
    # Arguments:
    
    # FUNCTION:
    
    # Settings:
    x = q
    xi = mu
    ans = NA
    
    # Univariate Case:
    if (is.vector(x) & dim == 1) {
        ans = pst(x, location = xi[1], scale = as.vector(Omega)[1], 
            shape = alpha[1], df = df)
    }
    
    # Multivariate Case:
    if (is.matrix(x)) {
        if (dim == ncol(x)) {
            ans = NULL
            for (i in 1:nrow(x) ) {
                ans = c(ans, pmst(x = x[i,], xi = xi, Omega = Omega, 
                    alpha = alpha, df = df))
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


rmvst = 
function(n, dim = 2, mu = rep(0, dim), Omega = diag(dim), 
alpha = rep(0, dim), df = 4)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Multivariate Skew Sudent-t Random Number Generator
    
    # Arguments:
    
    # FUNCTION:
    
    # Settings:
    ans = NA
    xi = mu
    
    # Univariate Case:
    if (dim == 1) {
        ans = as.matrix(rst(n, location = xi[1], 
            scale = as.vector(Omega)[1], shape = alpha[1], df = df))
    }
    
    # Multivariate Case:
    if (dim > 1) {
        ans = rmst(n, xi = xi, Omega = Omega, alpha = alpha, df = df)
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

