
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
# FUNCTION:                  DESCRIPTION:
#  BesselI                    Modified Bessel Function of first kind
#  BesselK                    Modified Bessel Function of third kind
#  BesselDI                   Derivative of BesselI
#  BesselDK                   Derivative of BesselK
# INTERNAL FUNCTION:         DESCRIPTION:
#  .BesselN                   For internal use only 
#  .Bessel01                   ...
#  .Bessel.MSTA1               ...
#  .Bessel.MSTA2               ...
#  .Bessel.ENVJ                ...
################################################################################


test.Bessel = 
function()
{
    # BesselI - Modified Bessel Function of first kind
    # BesselK - Modified Bessel Function of third kind

    # Modified Bessel Functions I and K, Abramowitz-Stegun, Chapter 9.6, p. 374
    
    # Abramowitz-Stegun, Figure 9.7, p. 374
    plot(x = c(0, 3), y = c(0, 2.5), type = "n", xlab = "x", 
        ylab = " I0  I1  K0  K1 ")
    grid()
    x = seq(0, 3, length = 301)
    lines(x, BesselI(x, 0))
    lines(x, BesselI(x, 1), lty = 3)
    lines(x, BesselK(x, 0))
    lines(x, BesselK(x, 1), lty = 3)
    
    # Abramowitz-Stegun, Figure 9.8, p. 375
    plot(x = c(0, 10), y = c(0, 1.9), type = "n", xlab = "x", ylab = "y")
    grid()
    x = seq(0, 10, length = 501)
    lines(x, exp(-x)*BesselI(x, 0))
    lines(x, exp(-x)*BesselI(x, 1), lty = 3)
    lines(x, exp(x)*BesselK(x, 0))
    lines(x, exp(x)*BesselK(x, 1), lty = 3)
    
    # Abramowitz-Stegun, Figure 9.9, p. 375
    # Use R's internally implemented Bessel Functions ...
    plot(x = c(-10, 10), y = c(-5, 30), type = "n", xlab = "nu", ylab = "y")
    grid()
    nu = seq(-10, 10, length = 501)
    lines(nu, besselI(x = 5, nu))
    lines(nu, besselK(x = 5, nu), lty = 3)

    # Bessel I0 and K0
    
    # Abramowitz-Stegun: Table 9.8, p. 416-422
    x = c(0.0, 0.01, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)  
    data.frame(x, I = exp(-x)*BesselI(x, 0), K = exp(x)*BesselK(x, 0)) 
    # Compare with R's internal function:
    data.frame(x, I = exp(-x)*besselI(x, 0), K = exp(x)*besselK(x, 0)) 
  
    # Bessel I2 and K2
    
    # Abramowitz-Stegun: Table 9.8, p. 416-422
    x = c(0.01, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)  
    data.frame(x, I = BesselI(x, 2)/x^2, K = BesselK(x, 2)*x^2) 
    # Compare with R's internal function:
    data.frame(x, I = besselI(x, 2)/x^2, K = besselK(x, 2)*x^2)
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.BesselAssociatedSeries = 
function()
{   
    # Associated Hyperbolic Series, Abramowitz-Stegun 9.6.39, p. 376
    besselI(x = 1, 0) + 2 * sum((-1)^(1:10)*besselI(x = 1, 2*(1:10)))
    1
    
    # Associated Hyperbolic Series, Abramowitz-Stegun 9.6.39, p. 376
    besselI(x = 1, 0) + 2 * sum(besselI(x = 1, 1:20))
    exp(1)
    
    # Associated Hyperbolic Series, Abramowitz-Stegun 9.6.39, p. 376
    besselI(x = 1, 0) + 2 * sum((-1)^(1:20)*besselI(x = 1, 1:20))
    exp(-1)  
    
    # Associated Hyperbolic Series, Abramowitz-Stegun 9.6.39, p. 376
    besselI(x = 1, 0) + 2 * sum(besselI(x = 1, 2*(1:10)))
    cosh(1)
    
    # Associated Hyperbolic Series, Abramowitz-Stegun 9.6.40, p. 376
    2* sum(besselI(x = 1, 2*(0:10)+1))
    sinh(1)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.BesselD = 
function()
{
    # BesselDI - Derivative of BesselI
    # BesselDK - Derivative of BesselK

    # Check:
    # I0'(x) = I1(x)
    x = c(0.01, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)
    BesselDI(x, 0) / BesselI(x, 1)
    
    # Check:
    # K0'(x) = -K1(x)
    x = c(0.01, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)
    - BesselDK(x, 0) / BesselK(x, 1)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


if (FALSE) {
    require(RUnit)
    testResult <- runTestFile("C:/Rmetrics/SVN/trunk/fOptions/tests/runit3D.R",
        rngKind = "Marsaglia-Multicarry", rngNormalKind = "Inversion")
    printTextProtocol(testResult)
}


################################################################################

