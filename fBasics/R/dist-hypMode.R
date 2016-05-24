
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
# FUNCTION:             DESCRIPTION:
#  hypMode               Computes the hyperbolic mode
# FUNCTION:             DESCRIPTION:
#  .hyp[1234]Mode         Internal functions called by 'hypMode'
################################################################################


hypMode <- 
function(alpha = 1, beta = 0, delta = 1, mu = 0, pm = c(1, 2, 3, 4))
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes the mode of the Hyperbolic PDF
  
    # FUNCTION:
    
    # Settings:
    pm = pm[1]
    
    # Return Value:
    ans = NA
    if (pm == 1) return(.hyp1Mode(alpha, beta, delta, mu))
    if (pm == 2) return(.hyp2Mode(alpha, beta, delta, mu))
    if (pm == 3) return(.hyp3Mode(alpha, beta, delta, mu))
    if (pm == 4) return(.hyp4Mode(alpha, beta, delta, mu))  
}


# ------------------------------------------------------------------------------


.hyp1Mode <- 
function(alpha = 1, beta = 0, delta = 1, mu = 0)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes the mode of the Hyperbolic PDF

    # FUNCTION:
    
    # Mode:
    ans = mu + delta * beta / sqrt(alpha^2 - beta^2)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.hyp2Mode <- 
function(zeta = 1, rho = 0, delta = 1, mu = 0)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes the hyperbolic mode in the 2nd parameterization

    # FUNCTION:
    
    # Parameter Change:
    alpha = zeta / ( delta * sqrt(1 - rho*rho) )
    beta = alpha * rho
    
    # Return Value:
    ans = hypMode(alpha, beta, delta, mu)
    ans
}


# ------------------------------------------------------------------------------


.hyp3Mode <-  
function(xi = 1/sqrt(2), chi = 0, delta = 1, mu = 0)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes the hyperbolic mode in the 3rd parameterization

    # FUNCTION:
    
    # Parameter Change:
    rho = chi / xi
    zeta = 1/xi^2 - 1   
    alpha = zeta / ( delta * sqrt(1 - rho*rho) )
    beta = alpha * rho
    ans = hypMode(alpha, beta, delta, mu)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.hyp4Mode <-  
function(a.bar = 1, b.bar = 0, delta  = 1, mu = 0)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes the hyperbolic mode in the 4th parameterization

    # FUNCTION:
    
    # Parameter Change:
    alpha = a.bar / delta
    beta = b.bar / delta
    ans = hypMode(alpha, beta, delta, mu)
    
    # Return Value:
    ans
}   


################################################################################

