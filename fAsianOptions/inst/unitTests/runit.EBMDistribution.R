
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
# FUNCTION:           EBM DENSITY APPROXIMATIONS:
#  dlognorm            log-Normal density an derivatives
#  plognorm            log-Normal, synonyme for plnorm
#  dgam                Gamma density, synonyme for dgamma
#  pgam                Gamma probability, synonyme for pgamma
#  drgam               Reciprocal-Gamma density
#  prgam               Reciprocal-Gamma probability
#  djohnson            Johnson Type I density
#  pjohnson            Johnson Type I probability
# FUNCTION :          MOMENTS FOR EBM DENSITY APPROXIMATIONS:
#  mnorm               Moments of Normal density
#  mlognorm            Moments of log-Normal density
#  mrgam               Moments of reciprocal-Gamma density
#  masian              Moments of Asian Option density
#  .DufresneMoments     Internal Function used by masian()
# FUNCTION:           NUMERICAL DERIVATIVES:
#  derivative          First and second numerical derivative
# FUNCTION:           ASIAN DENSITY:
#  d2EBM               Double Integrated EBM density
#  .thetaEBM            Internal Function used to compute *2EBM()
#  .psiEBM              Internal Function used to compute *2EBM()
#  dEBM                Exponential Brownian motion density
#  pEBM                Exponential Brownian motion probability              
#  .gxuEBM              Internal Function used to compute *EBM()
#  .gxtEBM              Internal Function used to compute *EBM()
#  .gxtuEBM             Internal Function used to compute *EBM()
#  dasymEBM            Exponential Brownian motion asymptotic density
################################################################################


test.lognorm = 
function()
{
    # dlognorm - log-Normal density an derivatives
    # plognorm - log-Normal, synonyme for plnorm

    # Calculate Log-Normal Density and its Derivatives:
    x = exp(seq(-2.8, 1.2, length = 100))
    y0 = dlognorm(x, deriv = 0)
    y1 = dlognorm(x, deriv = 1)
    y2 = dlognorm(x, deriv = 2) 
    
    # Compare with Numerical Differentiation:
    par(mfrow = c(2, 2))
    xa = exp(seq(-2.5, 1.5, length = 20))
    plot(x, y0, type = "l", main = "Log-Normal Density")
    plot(x, y1, type = "l", main = "1st Derivative")
    z = derivative(xa, dlognorm(xa, deriv = 0), deriv = 1)
    points(z$x, z$y, col = "steelblue")
    plot(x, y2, type = "l", main = "2nd Derivative")
    z = derivative(xa, dlognorm(xa, deriv = 0), deriv = 2)
    points(z$x, z$y, col = "steelblue")
    
    # Return Value:
    return()    
}
   

# ------------------------------------------------------------------------------


test.gam =
function()
{
    #  dgam - Gamma density, synonyme for dgamma
    #  pgam - Gamma probability, synonyme for pgamma

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.rgam = 
function()
{    
    # drgam - Reciprocal-Gamma density
    # prgam - Reciprocal-Gamma probability

    # Calculate Reciprocal-Gamma Density and its Derivaties:
    alpha = 2; beta = 1
    x = exp(seq(-2.8, 1.2, length = 100))
    y0 = drgam(x, alpha, beta, deriv = 0)
    y1 = drgam(x, alpha, beta, deriv = 1)
    y2 = drgam(x, alpha, beta, deriv = 2)
    
    # Compare with Numerical Differentiation:
    par(mfrow = c(2, 2))
    xa = exp(seq(-2.5, 1.5, length = 20))
    plot(x, y0, type = "l", main = "Rec-Gamma Density")
    plot(x, y1, type = "l", main = "1st Derivative")
    z = derivative(xa, drgam(xa, alpha, beta, deriv = 0), deriv = 1)
    points(z$x, z$y, col = "steelblue")
    plot(x, y2, type = "l", main = "2nd Derivative")
    z = derivative(xa, drgam(xa, alpha, beta, deriv = 0), deriv = 2)
    points(z$x, z$y, col = "steelblue")
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.johnson = 
function()
{        
    # djohnson - Johnson Type I density
    # pjohnson - Johnson Type I probability

    # Calculate Johnson-Type-I Density and its Derivaties:
    a = 0.3; b = 1.2; c = -0.2; d = 0.8
    x = exp(seq(-2.8, 1.2, length = 100))
    y0 = djohnson(x, a, b, c, d, deriv = 0)
    y1 = djohnson(x, a, b, c, d, deriv = 1)
    y2 = djohnson(x, a, b, c, d, deriv = 2)
    
    # Compare with Numerical Differentiation:
    par(mfrow = c(2, 2))
    xa = exp(seq(-2.5, 1.5, length = 20))
    plot(x, y0, type = "l", main = "Johnson Type I Density")
    plot(x, y1, type = "l", main = "1st Derivative")
    z = derivative(xa, djohnson(xa, a, b, c, d, deriv = 0), deriv = 1)
    points(z$x, z$y, col = "steelblue")
    plot(x, y2, type = "l", main = "2nd Derivative")
    z = derivative(xa, djohnson(xa, a, b, c, d, deriv = 0), deriv = 2)
    points(z$x, z$y, col = "steelblue")

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.moments = 
function()
{  
    # mnorm - Moments of Normal density
    # mlognorm - Moments of log-Normal density
    # mrgam - Moments of reciprocal-Gamma density
    # masian - Moments of Asian Option density
    # .DufresneMoments - Internal Function used by masian()

    # mnorm(mean = 0, sd = 1)
    mnorm()
    
    # mlognorm(meanlog = 0, sdlog = 1)
    mlognorm()
    
    # mrgam(alpha = 1/2, beta = 1)
    mrgam()                                                              # CHECK                                   
    
    # mjohnson(a, b, c, d)
    a = 0.3; b = 1.2; c = -0.2; d = 0.8
    mjohnson(a, b, c, d)                                                 # CHECK
    
    # masian(Time = 1, r = 0.045, sigma = 0.3) 
    masian()
    
    # .DufresneMoments(M = 4, Time = 1, r = 0.045, sigma = 0.30)
    .DufresneMoments(M = 12, Time = 1, r = 0.045, sigma = 0.30)
    
    # Return Value:
    return()    
}

    
# ------------------------------------------------------------------------------


test.d2EBM = 
function()
{  
    # d2EBM - Double Integrated EBM density
    # .thetaEBM - Internal Function used to compute *2EBM()
    # .psiEBM - Internal Function used to compute *2EBM()

    #  d2EBM(u, t = 1) 
    x = c(0.1, 0.5, 1, 2)
    
    # Density:
    d2 = d2EBM(u = x)
    d2
     
    # Compare with:
    d = dEBM(u = x)
    d
    
    # Print
    cbind(d2, d, difference = abs(d2-d))
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.dEBM = 
function()
{  
    #  dEBM - Exponential Brownian motion density
    #  pEBM - Exponential Brownian motion probability              
    #  .gxuEBM - Internal Function used to compute *EBM()
    #  .gxtEBM - Internal Function used to compute *EBM()
    #  .gxtuEBM - Internal Function used to compute *EBM()
    
    # Density:
    x = c(
        seq(-1.0, 0.0, length = 5), 
        seq( 0.0, 0.5, length = 20), 
        seq( 0.5, 5.0, length = 20))
    x = unique(sort(x))
    d = dEBM(u = x)
    print(d)
    par(mfrow = c(1,1))
    plot(x, y = d, type = "b", pch = 19, cex = 0.7)
    
    # Probability:
    x = c(-1, -0.5, 0, 0.1, 0.2, 0.5, 0.75, seq(1, 5, by = 0.5))
    p = pEBM(u = x)                          
    print(p)
    par(mfrow = c(1,1))
    plot(x, y = p, type = "b", pch = 19, cex = 0.7)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.dasymEBM = 
function()
{  
    # dasymEBM - Exponential Brownian motion asymptotic density
     
    # Density:
    x = c(
        seq(0.50, 1.10, length = 21),
        seq(1.10, 1.40, length = 21),
        seq(1.40, 5.00, length = 31))
    d = dasymEBM(u = x)
    print(d)
    par(mfrow = c(1,1))
    plot(x, y = d, type = "b", pch = 19, cex = 0.7)
    abline(h =0, lty = 3, col = "grey")
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


if (FALSE) {
    require(RUnit)
    testResult <- runTestFile("C:/Rmetrics/SVN/trunk/fOptions/tests/runit3A.R",
        rngKind = "Marsaglia-Multicarry", rngNormalKind = "Inversion")
    printTextProtocol(testResult)
}


################################################################################
