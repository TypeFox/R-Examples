
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
#  ghtMean                       Returns true GH Student-t mean
#  ghtVar                        Returns true GH Student-t variance
#  ghtSkew                       Returns true GH Student-t skewness
#  ghtKurt                       Returns true GH Student-t kurtosis
# FUNCTION:                     DESCRIPTION:
#  ghtMoments                    Returns true GH Student-t moments
################################################################################


ghtMean <-
function(beta=0.1, delta=1, mu=0, nu=10)
{
    # A function implemented by Diethelm Wuertz
    
    # Descriptions:
    #   Returns true Mean of the GH Student-t distribution
    
    # FUNCTION:
    
    # GH Parameters:
    alpha = abs(beta) + 1e-6
    lambda = -nu/2
    
    # Return Value:
    ghMean(alpha, beta, delta, mu, lambda)
}


# ------------------------------------------------------------------------------


ghtVar <- 
function(beta=0.1, delta=1, mu=0, nu=10)
{
    # A function implemented by Diethelm Wuertz
    
    # Descriptions:
    #   Returns true variance of the GH Student-t distribution
    
    # FUNCTION:
    
    # GH Parameters:
    alpha = abs(beta) + 1e-6
    lambda = -nu/2
    
    # Return Value:
    ghVar(alpha, beta, delta, mu, lambda)
}


# ------------------------------------------------------------------------------


ghtSkew <- 
function(beta=0.1, delta=1, mu=0, nu=10)
{
    # A function implemented by Diethelm Wuertz
    
    # Descriptions:
    #   Returns true Skewness of the GH Student-t distribution
    
    # FUNCTION:
    
    # GH Parameters:
    alpha = abs(beta) + 1e-6
    lambda = -nu/2
    
    # Return Value:
    ghSkew(alpha, beta, delta, mu, lambda)          
}


# ------------------------------------------------------------------------------


ghtKurt <- 
function(beta=0.1, delta=1, mu=0, nu=10)
{
    # A function implemented by Diethelm Wuertz
    
    # Descriptions:
    #   Returns true Kurtosis of the GH Student-t distribution
    
    # FUNCTION:
    
    # GH Parameters:
    alpha = abs(beta) + 1e-6
    lambda = -nu/2
    
    # Return Value:
    ghKurt(alpha, beta, delta, mu, lambda)
}


################################################################################


ghtMoments <-
function(order, type = c("raw", "central", "mu"),
    beta=0.1, delta=1, mu=0, nu=10)
{
    # A function implemented by Diethelm Wuertz
    
    # Descriptions:
    #   Returns moments of the GH Student-t distribution
    
    # FUNCTION:
    
    # Settings:
    type = match.arg(type)
    
    # GH Parameters:
    alpha = abs(beta) + 1e-6
    lambda = -nu/2
    
    # Return Value:
    ghMoments(order, type, alpha, beta, delta, mu, lambda)   
}


################################################################################

