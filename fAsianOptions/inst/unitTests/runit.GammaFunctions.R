
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
# FUNCTION:       DESCRIPTION:
#  erf             Error function
#  [gamma]         Gamma function 
#  [lgamma]        LogGamma function, returns log(gamma)
#  [digamma]       First derivative of of LogGamma, dlog(gamma(x))/dx
#  [trigamma]      Second derivative of of LogGamma, dlog(gamma(x))/dx
#  {tetragamma}    Third derivative of of LogGamma, dlog(gamma(x))/dx
#  {pentagamma}    Fourth derivative of LogGamma, dlog(gamma(x))/dx
#  [beta]*         Beta function
#  [lbeta]*        LogBeta function, returns log(Beta)
#  Psi             Psi(x) (Digamma) function
#  igamma          P(a,x) Incomplete Gamma Function
#  cgamma          Gamma function for complex arguments
#  Pochhammer      Pochhammer symbol
# NOTES:
#  Functions in [] paranthesis are part of the R's and SPlus' base distribution
#  Functions in {} paranthesis are only availalble in R
#  Function marked by []* are compute through the gamma function in SPlus
################################################################################


test.erf = 
function()
{
    # Error Function
    
    # erf(x) - Abramowitz-Stegun p. 310 ff
    erf(0.0)  
                     
    erf(0.5) - 0.5204998788   
    erf(1.0) - 0.8427007929
    erf(2.0) - 0.9953222650
    
    erf(-0.5) + 0.5204998788   
    erf(-1.0) + 0.8427007929
    erf(-2.0) + 0.9953222650
    
    x = seq(-5, 5, length = 101)
    y = erf(x)
    par(mfrow = c(1,1))
    plot(x, y, type = "b", pch = 19)

    # Symmetry Relation, p. 297
    # erf(-x) = - erf(x)
    erf(-0.5) + erf(0.5)
    erf(-pi)  + erf(pi)
    
    # Abramowitz-Stegun: Example 1, p. 304
    erf(0.745) - 0.707928920
  
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.erfc = 
function()
{
    # Complementary Error Function
    
    # Add functions:
    erfc = function(x) { 1 - erf(x) }
    erfc(0.5)
    erf(0.5) + erfc(0.5)
    
    # Abramowitz-Stegun: Example 2, p. 304
    1.1352e-11
    1 - erf(4.8)
  
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.gamma = 
function()
{
    # Gamma Function
    
    # Abramowitz-Stegun: Figure 6.1, p. 255
    x = seq(-4.01, 4.01, by = 0.011)
    
    # Plot:
    plot(x, gamma(x), ylim = c(-5,5), type = "l", main = "Gamma Function")
    grid()
    abline(h = 0, col = "red", lty = 3)
    for (i in c(-4:0)) abline(v = i, col = "white")
    for (i in c(-4:0)) abline(v = i, col = "red", lty = 3)
    
    # Add 1/Gamma:
    lines(x, 1/gamma(x), lty = 3)
  
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.Psi = 
function()
{    
    # Psi Function
    
    # Abramowitz-Stegun: Figure 6.2, p. 258
    x = seq(-4.01, 4.01, by = 0.011)
    plot(x, Psi(x), ylim = c(-5, 5), type = "l", main = "Psi Function")
    grid()
    abline(h = 0, col = "red", lty = 3)
    for (i in c(-4:0)) abline(v = i, col = "white")
    for (i in c(-4:0)) abline(v = i, col = "red", lty = 3)
  
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.incompleteGamma = 
function()
{        
    # Incomplete Gamma Function
    
    # Abramowitz-Stegun: Figure 6.3. 
    gammaStar = function(x, a) { igamma(x,a)/x^a }
    # ... create Figure as an exercise.

    # Abramowitz-Stegun: Formula 6.5.12
    # Relation to Confluent Hypergeometric Functions
    a = sqrt(2)
    x = pi
    Re ( (x^a/a) * kummerM(-x, a, 1+a) )
    Re ( (x^a*exp(-x)/a) * kummerM(x, 1, 1+a) )
    pgamma(x, a) * gamma(a)
    igamma(x, a)                                         
    
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.complexGamma = 
function()
{      
    # Complex Gamma Function
    
    # Abramowitz-Stegun: Tables 6.7, p. 277
    x = 1
    y = seq(0, 5, by = 1); x = rep(x, length = length(y))
    z = complex(real = x, imag = y)
    c = cgamma(z, log = TRUE)
    cbind(x, y, "Re ln" = Re(c), "Im ln" = Im(c))
    
    # Abramowitz-Stegun: Tables 6.7, p.287
    x = 2
    y = seq(0, 5, by = 1); x = rep(x, length = length(y))
    z = complex(real = x, imag = y)
    c = cgamma(z, log = TRUE)
    cbind(x, y, "Re ln" = Re(c), "Im ln" = Im(c))
    
    # cgamma -
    # Abramowitz-Stegun: Examples, p. 263
    options(digits = 10)
    # Example 1, 2
    gamma(6.38); lgamma(56.38)                            
    # 3, 4
    Psi(6.38); Psi(56.38)                                 
    # 5
    cgamma(complex(real = 1, imag = -1), log = TRUE )     
    # 6
    cgamma(complex(real = 1/2, imag = 1/2), log = TRUE )  
    # 7, 8
    cgamma(complex(real = 3, imag = 7), log = TRUE )       

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.Pochhammer = 
function()
{      
    # Pochhammer Symbol
    
    # Abramowitz-Stegun: Formula 6.1.22, p. 256
    
    # Pochhammer(x, n)
    Pochhammer(x = 1, n = 0) - 1
    Pochhammer(x = 1, n = 1) - 1
    Pochhammer(x = 2, n = 2) - gamma(2+2)/gamma(2)
        
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------

   
if (FALSE) {
    require(RUnit)
    testResult <- runTestFile("C:/Rmetrics/SVN/trunk/fOptions/tests/runit3B.R",
        rngKind = "Marsaglia-Multicarry", rngNormalKind = "Inversion")
    printTextProtocol(testResult)
}


################################################################################

