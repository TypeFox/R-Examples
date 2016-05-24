
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
# FUNCTION:                  ELLIPTICAL COPULAE PARAMETER FITTING:
#  ellipticalCopulaSim        Simulates bivariate elliptical copula
#  ellipticalCopulaFit        Fits the paramter of an elliptical copula
################################################################################


test.copulaSim = 
function()
{  
    # Arguments:
    # ellipticalCopulaSim(n, rho = 0.75, param = NULL, 
    #   type = c("norm", "cauchy", "t")) 

    # Normal Copula:
    rho = 0.5
    R = ellipticalCopulaSim(n = 1000, rho = rho)
    R[1:10, ]
    plot(R, pch = 19)

    # Cauchy Copula:
    rho = runif(1, -1, 1)
    R = ellipticalCopulaSim(n = 100, rho = rho, type = "cauchy")
    R[1:10, ]
    plot(R, pch = 19)
    
    # Student-t Copula:
    rho = runif(1, -1, 1)
    nu = runif(1, 3, 20)
    print(c(rho, nu))
    R = ellipticalCopulaSim(n = 1000, rho = rho, param = nu, type = "t")
    R[1:10, ]
    plot(R, pch = 19)
    
    # The remaining Copulae are not yet implemented ...
    
    # Return Value:
    return()
}
    
    
# ------------------------------------------------------------------------------


test.copulaFit = 
function()
{     
    # Arguments:
    # ellipticalCopulaFit(u, v = NULL, type = c("norm", "cauchy", "t"), ...) 
  
    # Fit Normal Copula:
    rho = 0.5
    R = ellipticalCopulaSim(n = 1000, rho = rho)
    fit = ellipticalCopulaFit(u = R[,1], v = R[,2])
    fit
    rho - fit$par
    plot(c(-1,1), c(-1,1), xlab = "rho", ylab = "estimate", main = "Normal")
    for ( i in 1:100) {
        rho = runif(1, -1, 1)
        R = ellipticalCopulaSim(n = 1000, rho = rho)
        fit = ellipticalCopulaFit(R)
        points(rho, fit$par)
        print(c(rho, fit$par))
    }

    # Fit Cauchy Copula:
    rho = runif(1, -1, 1)
    R = ellipticalCopulaSim(n = 100, rho = rho, type = "cauchy")
    ellipticalCopulaFit(R, type = "cauchy")
    rho
    plot(c(-1,1), c(-1,1), main = "Cauchy")
    for ( i in 1:100) {
        rho = runif(1, -1, 1)
        R = ellipticalCopulaSim(n = 1000, rho = rho, type = "cauchy")
        fit = ellipticalCopulaFit(R, type = "cauchy")
        points(rho, fit$par)
        print(c(rho, fit$par))
    }
    
    # Fit Student-t Copula:
    rho = runif(1, -1, 1)
    nu = runif(1, 3, 20)
    print(c(rho, nu))
    R = ellipticalCopulaSim(n = 1000, rho = rho, param = nu, type = "t")
    ellipticalCopulaFit(R, type = "t")
    plot(c(-1,1), c(-1,1), main = "Student-t")
    for ( i in 1:100) {
        rho = runif(1, -1, 1)
        nu = runif(1, 3, 20)
        R = ellipticalCopulaSim(n = 1000, rho = rho, param = nu, type = "t")
        fit = ellipticalCopulaFit(R, type = "t")
        points(rho, fit$par[1])
        print(c(rho, nu, fit$par))
    }
    
    # The remaining Copulae are not yet implemented ...
    
    # Return Value:
    return()
}
        
  
################################################################################
   
