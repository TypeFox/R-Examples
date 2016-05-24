
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
# You should have received A copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:            DESCRIPTION:
#  dsnig                Returns density of the SNIG distribution
#  psnig                Returns probabilities of the SNIG distribution
#  qsnig                Returns quantiles of the SNIG distribution
#  rsnig                Generates SNIG distributed random variates
# FUNCTION:            DESCRIPTION:
#  .psnigC              Fast psnig from C code
#  .qsnigC              Fast qsnig from C code
################################################################################


dsnig <-  
function(x, zeta = 1, rho = 0, log = FALSE) 
{
    # Description:
    #   Returns density of the snig distribution
    
    # FUNCTION:
 
    # Parameters:
    if (length(zeta) == 2) {
       rho = zeta[2]
       zeta = zeta[1]
    } 
    
    # Compute Density - Quick and Dirty:
    ans = dsgh(x, zeta, rho, lambda = -0.5, log = log)
    
    # Log:
    if(log) ans = log(ans)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


psnig <-  
function(q, zeta = 1, rho = 0) 
{
    # Description:
    #   Returns probabilities of the snig distribution
    
    # FUNCTION:
    
    # Compute Probabilities - Quick and Dirty:
    ans = psgh(q, zeta, rho, lambda = -0.5)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


qsnig <-  
function(p, zeta = 1, rho = 0) 
{
    # Description:
    #   Returns quantiles of the snig distribution
    
    # FUNCTION:
    
    # Compute Quantiles:
    ans = qsgh(p, zeta, rho, lambda = -0.5)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


rsnig <-  
function(n, zeta = 1, rho = 0) 
{
    # Description:
    #   Generates snig distributed random variates
    
    # FUNCTION:
    
    # Generate Random Numbers:
    ans = rsgh(n, zeta, rho, lambda = -0.5)
    
    # Return Value:
    ans
}


################################################################################


.psnigC <-  
function(q, zeta = 1, rho = 0) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns quantiles of the snig distribution
    
    # FUNCTION:
      
    # Compute Quantiles:
    param = .paramGH(zeta, rho, lambda = -0.5)
    ans = .pnigC(q, param[1], param[2], param[3], param[4])
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.qsnigC <-  
function(p, zeta = 1, rho = 0) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns quantiles of the snig distribution
    
    # FUNCTION:
      
    # Compute Quantiles:
    param = .paramGH(zeta, rho, lambda = -0.5)
    ans = .qnigC(p, param[1], param[2], param[3], param[4])
    
    # Return Value:
    ans
}


################################################################################

