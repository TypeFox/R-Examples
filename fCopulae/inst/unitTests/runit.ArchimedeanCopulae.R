
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
# FUNCTION:                  ARCHIMEDEAN COPULAE RANDOM VARIATES:
#  rarchmCopula               Generates Archimedean copula random variates
#  rarchmSlider               Displays interactively archimedean probability
# FUNCTION:                  ARCHIMEDEAN COPULAE PROBABILITY:
#  parchmCopula               Computes Archimedean copula probability 
#  parchmSlider               Displays interactively archimedean probability 
# FUNCTION:                  ARCHIMEDEAN COPULAE DENSITY:
#  darchmCopula               Computes Archimedean copula density 
#  darchmSlider                Displays interactively archimedean density 
# FUNCTION:                  SPECIAL BIVARIATE COPULA:
#  rgumbelCopula              Generates fast gumbel random variates
#  pgumbelCopula              Computes bivariate Gumbel copula probability
#  dgumbelCopula              Computes bivariate Gumbel copula density
################################################################################

    
test.rarchmCopula = 
function()
{
    # Arguments:
    # rarchmCopula(n, alpha = NULL, type = archmList())
    
    # Random Variates - Check all Types:
    for (type in archmList()) {
        R = rarchmCopula(n = 5, alpha = NULL, type = type)
        cat("\n")
        print(type)
        print(R)
    }
      
    # Return Value:
    return()    
} 


# ------------------------------------------------------------------------------  
    

test.rarchmSlider = 
function()
{        
    # Arguments:
    # rarchmSlider(B = 10)
    
    # Try Slider:
    # rarchmSlider()
    NA
  
    # Return Value:
    return()    
}  


################################################################################  
    

test.parchmCopula = 
function()
{
    # Arguments:
    # parchmCopula(u = 0.5, v = u, alpha = NULL, type = archmList(), 
    #   output = c("vector", "list"), alternative = FALSE)
    
    # u - single input value:
    parchmCopula()
    parchmCopula(0.5)
    parchmCopula(0.5, 0.25)
    
    # u - input vector:
    U = (0:10)/10
    V = U
    parchmCopula(U)
    parchmCopula(u = U, v = V)
    parchmCopula(u = U, v = rev(V))
    
    # u - input matrix:
    parchmCopula(cbind(U, V))
    
    # u - input list:
    u = grid2d()
    u
    parchmCopula(u) # output = "vector"
    parchmCopula(u, output = "list")
    diff = parchmCopula(u) - parchmCopula(u, alternative = TRUE)
    mean(abs(diff))
   
    # Check All Types:
    u = grid2d()
    for (type in paste(1:22)) {
        cop1 = parchmCopula(u, type = type, output = "list")
        cop2 = parchmCopula(u, type = type, output = "list", alternative = TRUE)
        cat("Type: ", type, "\t Difference: ", mean(abs(cop1$z-cop2$z)), "\n")
        persp(cop1, main = type, theta = -40, phi = 30, col = "steelblue")
    }  
  
    # Return Value:
    return()    
} 


# ------------------------------------------------------------------------------  
    

test.parchmSlider = 
function()
{        
    # Arguments:
    # parchmSlider(type = c("persp", "contour"), B = 10)
    
    # Try Perspective Slider:
    # parchmSlider()
    NA 
    
    # Try Contour Slider:
    # parchmSlider("contour")
    NA
    
    # Return Value:
    return()    
}  


################################################################################  
    

test.darchmCopula = 
function()
{
    # Arguments:
    # darchmCopula(u = 0.5, v = u, alpha = NULL, type = archmList(), 
    #   output = c("vector", "list"), alternative = FALSE)
    
    # u - single input value:
    darchmCopula()
    darchmCopula(0.5)
    darchmCopula(0.5, 0.25)
    
    # u - input vector:
    U = (0:10)/10
    V = U
    darchmCopula(U)
    darchmCopula(u = U, v = V)
    darchmCopula(u = U, v = rev(V))
    
    # u - input matrix:
    darchmCopula(cbind(U, V))
    
    # u - input list:
    u = grid2d()
    u
    darchmCopula(u) # output = "vector"
    darchmCopula(u, output = "list")
    
    # Check All Types:
    u = grid2d(x = (0:25)/25)
    for (type in archmList()) {
        cop1 = darchmCopula(u, type = type, output = "list")
        cop2 = darchmCopula(u, type = type, output = "list", 
            alternative = TRUE)
        diff = abs(cop1$z-cop2$z)
        diff = diff[!is.na(diff)]
        cat("Type: ", type, "\t Difference: ", mean(diff), "\n")
        persp(cop2, main = type, theta = -40, phi = 30, col = "steelblue")
    }
  
    # Return Value:
    return()    
} 

    
# ------------------------------------------------------------------------------  
    

test.darchmSlider = 
function()
{        
    # Arguments:
    # darchmSlider(type = c("persp", "contour"), B = 10)
    
    # Try Perspective Slider:
    # darchmSlider()
    NA
    
    # Try Contour Slider:
    # darchmSlider("contour")
    NA
    
    # Return Value:
    return()    
}  


################################################################################



test.rgumbelCopula = 
function()
{              
    # Generates fast gumbel random variates

    # Copula:
    rgumbelCopula()
    
    # Return Value:
    return()    
} 

    
# ------------------------------------------------------------------------------  
    

test.pgumbelCopula = 
function()
{              
    # Computes bivariate Gumbel copula probability

    # Copula:
    pgumbelCopula()
    
    # Return Value:
    return()    
} 

    
# ------------------------------------------------------------------------------  
    

test.dgumbelCopula = 
function()
{              
    # Computes bivariate Gumbel copula density
  
    # Copula:
    dgumbelCopula()
    
    # Return Value:
    return()    
} 

  
################################################################################

