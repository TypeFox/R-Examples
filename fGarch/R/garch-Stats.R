
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
# FUNCTION:                DESCRIPTION: 
#  garchKappa               Computes Expection for APARCH Models
#  .garchKappaFun           Internal function used by garchKappa()
# FUNCTION:                DESCRIPTION:
#  .truePersistence         Computes true persistence
################################################################################


garchKappa <-  
    function(cond.dist = c("norm", "ged", "std", "snorm", "sged", "sstd",
    "snig"), gamma = 0, delta = 2, skew = NA, shape = NA)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Expection for APARCH Models
    
    # FUNCTION:
    
    # Compute kappa:
    kappa = integrate(.garchKappaFun, lower = -Inf, upper = Inf, cond.dist = 
        cond.dist[1], gamma = gamma, delta = delta, skew = skew, shape = 
        shape)[[1]] 
    names(kappa) = "kappa"
    attr(kappa, "control") = 
        c(gamma = gamma, delta = delta, skew = skew, shape = shape)
    attr(kappa, "cond.dist") = cond.dist[1]
    
    # Return Value:
    kappa
}


# ------------------------------------------------------------------------------


.garchKappaFun <-  
    function(x, 
    cond.dist = c("norm", "ged", "std", "snorm", "sged", "sstd", "snig"), 
    gamma = 0, delta = 2, skew = NA, shape = NA)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Internal function used by kappa()

    # FUNCTION:
    
    # Compute Expectation Value for ...
    funcE = (abs(x) - gamma*x)^delta
    
    # Select Appropriate Conditional Density:
    cond.dist = cond.dist[1]
    if (cond.dist == "norm") {
        fun = funcE * dnorm(x)
    }
    if (cond.dist == "ged") {
        fun = funcE * dged(x, nu = shape) 
    }
    if (cond.dist == "std") {
        fun = funcE * dstd(x, nu = shape) 
    }
    if (cond.dist == "snorm") {
        fun = funcE * dsnorm(x, xi = skew)
    }
    if (cond.dist == "sged") {
        fun = funcE * dsged(x, nu = shape, xi = skew) 
    }
    if (cond.dist == "sstd") {
        fun = funcE * dsstd(x, nu = shape, xi = skew) 
    }
    if (cond.dist == "snig") {
        fun = funcE * dsnig(x, zeta = shape, rho = skew) 
    }
    
    # Return Value:
    fun
} 


################################################################################


.truePersistence <- 
    function(fun = "norm", alpha = 1, gamma = 0, beta = 0, delta = 1, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes persistence for an APARCH process
    
    # Arguments:
    #   fun - name of density functions of APARCH innovations
    #   alpha, gamma - numeric value or vector of APARCH coefficients,
    #       must be of same length  
    #   beta - numeric value or vector of APARCH coefficients
    #   delta - numeric value of APARCH exponent
    
    # Note:
    #   fun is one of: norm, snorn, std, sstd, ged, sged, snig
    
    # FUNCTION:
    
    # Match Density Function:
    fun = match.fun(fun)
    
    # Persisgtence Function: E(|z|-gamma z)^delta
    e = function(x, gamma, delta, ...) {
        (abs(x)-gamma*x)^delta * fun(x, ...)
    }
        
    # Compute Persistence by Integration:
    persistence = sum(beta)
    for (i in 1:length(alpha)) {
        I = integrate(e, -Inf, Inf, subdivisions = 1000, 
            rel.tol = .Machine$double.eps^0.5, 
            gamma = gamma[i], delta = delta, ...)
        persistence = persistence + alpha[i] * I[[1]]
    }
    
    # Warning:
    if (persistence >= 1) {  
        p = as.character(round(persistence, digits = 3))
        warning(paste("Divergent persistence p =", p))
    }
    
    # Return Value:
    persistence
}


################################################################################

