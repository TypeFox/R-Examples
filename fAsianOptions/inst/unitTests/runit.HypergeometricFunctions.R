
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
# FUNCTION:       KUMMER DESCRIPTION:               
#  kummerM         Computes Confluent Hypergeometric Function of the 1st Kind
#  kummerU         Computes Confluent Hypergeometric Function of the 2nd Kind
# FUNCTION:       WHITTAKER DESCRIPTION:
#  whittakerM      Computes Whittaker's M Function
#  whittakerW      Computes Whittaker's M Function
# FUNCTION:       HERMITE POLYNOMIAL:
#  hermiteH        Computes the Hermite Polynomial
################################################################################


test.kummer = 
function()
{
    # kummerM(x, a, b, lnchf = 0, ip = 0) 
    # kummerU(x, a, b, ip = 0)
    
    # Relation to Modified Bessel Function:
    # Abramowitz-Stegun: Formula 13.6.3, p. 509
    x = c(0.001, 0.01, 0.1, 1, 10, 100, 1000)  
    nu = 1; a = nu+1/2; b = 2*nu+1
    M = Re ( kummerM(x = 2*x, a = a, b = b) )
    Bessel = gamma(1+nu) * exp(x)*(x/2)^(-nu) * BesselI(x, nu)
    cbind(x, M, Bessel) 
    
    # Relation to Hyperbolic Function:
    # Abramowitz-Stegun: Formula 13.6.14, p. 509
    x = c(0.001, 0.01, 0.1, 1, 10, 100, 1000)  
    M = Re ( kummerM(2*x, a = 1, b = 2) )
    Sinh = exp(x)*sinh(x)/(x)
    cbind(x, M, Sinh)
    
    # Relation to Complex Hyperbolic Function:
    # Abramowitz-Stegun: Formula 13.6.14, p. 509
    y = rep(1, length = length(x))
    x = complex(real = x, imag = y)
    M = kummerM(2*x, a = 1, b = 2)
    Sinh = exp(x)*sinh(x)/(x)
    cbind(x, M, Sinh)
    
    # Abramowitz-Stegun: Examples, p. 511
    M = function(a, b, x) { Re (kummerM(x, a, b)) }     
    # Example 1
    M(  0.3,  0.2, -0.1)  - 0.8578490                       
    # 2
    M( 17.0, 16.0, 1)     - 2.8881744                       
    # 3
    M( -1.3,  1.2, 0.1)   - 0.8924108                      
    # 4
    # M( -1, -1.0, 0)     # undefined                   
    # 6
    M( 0.9,   0.1, 10)    - 1227235                       
    # 7
    M(-52.5,  0.1, 1)     - 16.34                                        # CHECK
    
    # Abramowitz-Stegun: Examples, p. 511
    U = function(a, b, x) { Re (kummerU(x, a, b)) }     
    # 9
    U( 1.1, 0.2, 1)     - 0.38664                           
    U(-0.9, 0.2, 1)     - 0.91272                         
    # 10
    U( 0.1, 0.2, 1)*0.9 - 0.85276                            
    # 11
    U( 1.0, 0.1, 100)   - 0.0098153                                      # CHECK
    # 12
    U( 0.1, 0.2, 0.01)  - 1.09                        
    
    # Abramowitz-Stegun: Example 17, Figure 13.2, p. 513
    # M(-4.5, 1, x)
    x = seq(0, 17, length = 200)
    plot(x = x, y = kummerM(x, -4.5, 1), type = "l", ylim = c(-25,125),
        main = "Figure 13.2:  M(-4.5, 1, x)")
    lines(x = c(0, 16), y = c(0, 0), col = 2)
    
    # Abramowitz-Stegun: Example 17, Figure 13.3, p. 513
    # M(a, 1, x)
    x = seq(0, 8, length = 200)
    plot(x = c(0, 8), y = c(-10, 10), type = "n", 
        main = "Figure 13.2:  M(-4.5, 1, x)")
    grid()
    for (a in seq(-4, 0, by = 0.5))
        lines(x = x, y = kummerM(x, a, 1), type = "l")
    abline(h = 0, lty = 3, col = "red")
    
    # Abramowitz-Stegun: Example 17, Figure 13.4, p. 513
    # M(a, 0.5, x)
    x = seq(0, 7, length = 200)
    plot(x = c(0, 7), y = c(-10, 15), type = "n", 
        main = "Figure 13.2:  M(-4.5, 1, x)")
    grid()
    for (a in seq(-4, 0, by = 0.5))
        lines(x = x, y = kummerM(x, a, 0.5), type = "l")
    abline(h = 0, lty = 3, col = "red")
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.whittaker = 
function()
{
    # whittakerM(x, kappa, mu, ip = 0) 
    # whittakerW(x, kappa, mu, ip = 0)
     
    # Abramowitz-Stegun: Example 13
    AS = c(1.10622, 0.57469)
    W = c(
     whittakerM(x = 1, kappa = 0, mu = -0.4),
     whittakerW(x = 1, kappa = 0, mu = -0.4) )
    data.frame(AS, W)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.hermite = 
function()
{
    # Hermite Polynomial - internally computed from Kummer U
    
    # hermiteH(x, n, ip = 0)

    # http://mathworld.wolfram.com/HermitePolynomial.html
    x = seq(-2, 2, length = 401)
    par(mfrow = c(1, 1))
    plot(x = c(-2,2), y = c(-30, 30), type = "n", main = "Hermite Polynomials")
    grid()
    for (i in 1:4) lines(x, hermiteH(x, i), col = i)
    
    # Test H4:
    H4 = function(x) { 16*x^4 -48*x^2 + 12 }
    x = -5:5
    cbind(x, hermite = hermiteH(x, 4), H = H4(x))
    
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


if (FALSE) {
    require(RUnit)
    testResult <- runTestFile("C:/Rmetrics/SVN/trunk/fOptions/tests/runit3C.R",
        rngKind = "Marsaglia-Multicarry", rngNormalKind = "Inversion")
    printTextProtocol(testResult)
}


################################################################################

