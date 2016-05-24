
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


################################################################################
# FUNCTION:                     DESCRIPTION:
#  nigMean                       Returns true NIG mean 
#  nigVar                        Returns true NIG variance 
#  nigSkew                       Returns true NIG skewness 
#  nigKurt                       Returns true NIG kurtosis 
# FUNCTION:                     DESCRIPTION:
#  nigMoments                    Returns true NIG moments 
# FUNCTION:                     DESCRIPTION:
#  nigShapeTriangle              Plots NIG Shape Triangle
################################################################################


nigMean <- 
function(alpha = 1, beta = 0, delta = 1, mu = 0)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns true NIG mean 
    
    # FUNCTION:
    
    # Mean:
    gamma = sqrt(alpha^2 - beta^2)
    ans = c(mean = (mu + delta * beta / gamma)[[1]])
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


nigVar <- 
function(alpha = 1, beta = 0, delta = 1, mu = 0)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns true NIG variance
    
    # FUNCTION:
    
    # Variance:
    gamma = sqrt(alpha^2 - beta^2)
    ans = c(var = (delta * alpha^2 / gamma^3)[[1]])
    
    # Return Value:
    ans
}



# ------------------------------------------------------------------------------


nigSkew <- 
function(alpha = 1, beta = 0, delta = 1, mu = 0)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns true NIG skewness
    
    # FUNCTION:
    
    # Skewness:
    gamma = sqrt(alpha^2 - beta^2)
    ans = c(skew = (3*beta / ( alpha * sqrt(delta*gamma) ))[[1]] )
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


nigKurt <- 
function(alpha = 1, beta = 0, delta = 1, mu = 0)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns true NIG kurtosis
    
    # FUNCTION:
    
    # Skewness:
    gamma = sqrt(alpha^2 - beta^2)
    ans = c(kurt = (3 * ( 1 + 4 * beta^2 / alpha^2) / (delta * gamma))[[1]])
    
    # Return Value:
    ans
}


###############################################################################


nigMoments <-
function(order, type = c("raw", "central", "mu"), 
    alpha=1, beta=0, delta=1, mu=0)
{
    # A function implemented by Diethelm Wuertz
    
    # Descriptions:
    #   Returns true moments of the NIG distribution
    
    # FUNCTION:
    
    # Settings:
    type = match.arg(type)
    
    # Moments:
    lambda = -1/2
    if (type == "raw") {
        ans = .ghRawMoments(order, alpha, beta, delta, mu, lambda)
    } else if (type == "central") {
        ans = .ghCentralMoments(order, alpha, beta, delta, mu, lambda)
    } else if (type == "central") {
        ans = .ghMuMoments(order, alpha, beta, delta, mu, lambda)  
    }
    
    # Return Value:
    ans   
}


################################################################################

   
nigShapeTriangle <- 
function(object, add = FALSE, labels = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots NIG Shape Triangle
    
    # Arguments:
    #   object - an object of class 'fDISTFIT' as returned by the
    #   function nigFit()
    
    # Example:
    #   nigShapeTriangle(nigFit(rnig(100), doplot = FALSE))
    
    # FUNCTION:
    
    # Settings:
    stopifnot(class(object) == "fDISTFIT")
    
    # Plot Frame:
    if (labels) {
        xlab = "Asymmetry: chi"
        ylab = "Steepness: zeta"
        main = "NIG Shape Traingle"
    } else {
        xlab = ylab = main = ""
    }
    if (!add) {
        x = c(-1, 0, 1, -1)
        y = c( 1, 0, 1,  1)
        plot(x, y, type = "l", xlab = xlab, ylab = ylab, main = main, ...)
        for (s in c(0.8, 0.6, 0.4, 0.2)) 
            lines(c(-s, s), c(s, s), lty = 3, col = "grey")  
        lines(c(0, 0), c(0, 1), lty = 3, col = "grey")
    }
    
    # Transform:
    par = object@fit$estimate
    names(par) = c("alpha", "beta", "delta", "mu")
    alpha = par[1] 
    beta = par[2]
    delta = par[3]
    mu = par[4]   
    
    # Add Points:
    zeta = 1/sqrt(1+delta*sqrt(alpha^2-beta^2))
    chi = zeta*(beta/alpha)
    if (labels) {
        points(chi, zeta, pch = 19, ...)
    } else {
        points(chi, zeta, ...)
    }   
    
    # Result:
    ans = list(chi = chi[[1]], zeta = zeta[[1]])
    attr(ans, "control")<-par
    
    # Return Value:
    ans
}


################################################################################

