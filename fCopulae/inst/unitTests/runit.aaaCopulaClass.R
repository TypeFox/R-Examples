
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
# FUNCTION:                 SPECIFICATION:
#  fCOPULA                   S4 class representation          
#  show                      S4 print method for copula specification
# FUNCTION:                 FRECHET COPULAE:
#  pfrechetCopula            Computes Frechet copula probability
# FUNCTION:                 SPEARMAN'S RHO:
#  .copulaRho                Spearman's rho by integration for "ANY" copula 
################################################################################


test.pfrechetCopula = 
function()
{
    # pfrechetCopula(u = 0.5, v = u, type = c("m", "pi", "w"), 
    #   output = c("vector", "list")) 
 
    # Grid:
    grid2d()
    grid2d(x = (0:10)/10)   
    
    # Vector - M Copula:
    copula.vector = pfrechetCopula(u = grid2d(), output = "vector", type = "m")
    copula.vector
    class(copula.vector)
    cbind(u = grid2d()$x, v = grid2d()$y, C = copula.vector)
    
    # List - M Copula:
    copula.list = pfrechetCopula(u = grid2d(), output = "list", type = "m")
    copula.list
    class(copula.list)
    persp(copula.list, theta = -40, phi = 30, col = "steelblue", ps = 9)
    
    # Vector - Pi Copula:
    copula.vector = pfrechetCopula(u = grid2d(), output = "vector", type = "pi")
    copula.vector
    class(copula.vector)
    cbind(u = grid2d()$x, v = grid2d()$y, C = copula.vector)
    
    # List - Pi Copula:
    copula.list = pfrechetCopula(u = grid2d(), output = "list", type = "pi")
    copula.list
    class(copula.list)
    persp(copula.list, theta = -40, phi = 30, col = "steelblue", ps = 9)
    
    # Vector - W Copula:
    copula.vector = pfrechetCopula(u = grid2d(), output = "vector", type = "w")
    copula.vector
    class(copula.vector)
    cbind(u = grid2d()$x, v = grid2d()$y, C = copula.vector)
    
    # List - W Copula:
    copula.list = pfrechetCopula(u = grid2d(), output = "list", type = "w")
    copula.list
    class(copula.list)
    persp(copula.list, theta = -40, phi = 30, col = "steelblue", ps = 9)
    
    # Return Value:
    return()    
}


################################################################################


test.copulaRho = 
function()
{
    # .copulaRho(rho = NULL, alpha = NULL, param = NULL, 
    #     family = c("elliptical", "archm", "ev", "archmax"), 
    #     type = NULL, error = 1e-3, ...)
    
    # Elliptical:
    fCopulae:::.copulaRho(rho = 0.5, family = "elliptical", type = "norm")
    
    # Return Value:
    return()    
}

  
################################################################################

  