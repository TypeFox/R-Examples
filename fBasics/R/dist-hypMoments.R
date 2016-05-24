
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
#  hypMean                       Returns true hyperbolic mean
#  hypVar                        Returns true hyperbolic variance
#  hypSkew                       Returns true hyperbolic skewness
#  hypKurt                       Returns true hyperbolic kurtosis
# FUNCTION:                     DESCRIPTION:
#  hypMoments                    Returns true hyperbolic moments
################################################################################


hypMean <-
function(alpha=1, beta=0, delta=1, mu=0)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns true hyperbolic mean
    
    # FUNCTION:
    
    # Return Value:
    ghMean(alpha, beta, delta, mu, lambda=1)
}


# ------------------------------------------------------------------------------


hypVar <- 
function(alpha=1, beta=0, delta=1, mu=0)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns true hyperbolic variance
    
    # FUNCTION:
    
    # Return Value:
    ghVar(alpha, beta, delta, mu, lambda=1)
}


# ------------------------------------------------------------------------------


hypSkew <- 
function(alpha=1, beta=0, delta=1, mu=0)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns true hyperbolic skewness
    
    # FUNCTION:
    
    # Return Value:
    ghSkew(alpha, beta, delta, mu, lambda=1)      
}


# ------------------------------------------------------------------------------


hypKurt <- 
function(alpha=1, beta=0, delta=1, mu=0)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns true hyperbolic kurtosis
    
    # FUNCTION:
   
    # Return Value:
    ghKurt(alpha, beta, delta, mu, lambda=1)
}


# ------------------------------------------------------------------------------


hypMoments <-
function(order, type = c("raw", "central", "mu"), 
    alpha=1, beta=0, delta=1, mu=0)
{
    # A function implemented by Diethelm Wuertz
    
    # Descriptions:
    #   Returns true moments of the hyperbolic distribution
    
    # FUNCTION:
    
    # Settings:
    type = match.arg(type)
    
    # Moments:
    lambda = 1
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

