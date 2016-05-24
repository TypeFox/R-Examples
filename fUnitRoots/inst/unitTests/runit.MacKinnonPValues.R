
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
# FUNCTION:          MC KINNON'S PROBABILIY AND QUANTILES:
#  punitroot          Returns cumulative probability for unit root distributions
#  qunitroot          Returns quantiles for unit root distributions
#  unitrootTable      Returns McKinnon's unitroot finite sample test table
# FUNCTION:          INTERNAL UTILITY FUNCTIONS:
#  .strsplitUrcval    Implements string split function for S-Plus compatibility
#  .urcval            Implements unit root statists
#  .probsUrcval       Implements probability values
# FUNCTION:          INTERNAL DATA SETS:
#  .urc1 ... .urc12   Statistical values for unitroot data
################################################################################


test.asymptoticUnitroot = 
function()
{
    # n.sample = 0 | Inf
    # trend = c("c", "nc", "ct", "ctt"), 
    # statistic = c("t", "n")

    # Asymptotic quantiles
    tol = .Machine$double.eps^0.25
    X = c(0.05, 0.10, 0.50, 0.90, 0.95)
    
    Q = qunitroot(X, trend = "c",   statistic = "t")
    P = punitroot(Q,    trend = "c",   statistic = "t")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)
    
    Q = qunitroot(X, trend = "nc",  statistic = "t")
    P = punitroot(Q,    trend = "nc",  statistic = "t")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)
    
    Q = qunitroot(X, trend = "ct",  statistic = "t")
    P = punitroot(Q,    trend = "ct",  statistic = "t")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)
    
    Q = qunitroot(X, trend = "ctt", statistic = "t")
    P = punitroot(Q,    trend = "ctt", statistic = "t")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol) 
    
    Q = qunitroot(X, trend = "c",   statistic = "n")
    P = punitroot(Q,    trend = "c",   statistic = "n")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)
    
    Q = qunitroot(X, trend = "nc",  statistic = "n")
    P = punitroot(Q,    trend = "nc",  statistic = "n")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)
    
    Q = qunitroot(X, trend = "ct",  statistic = "n")
    P = punitroot(Q,    trend = "ct",  statistic = "n")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)
    
    Q = qunitroot(X, trend = "ctt", statistic = "n")
    P = punitroot(Q,    trend = "ctt", statistic = "n")   
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.finiteSizeUnitroot = 
function()
{    
    # Finite size quantiles - Sample Size = 100
    tol = .Machine$double.eps^0.25
    X = c(0.05, 0.10, 0.50, 0.90, 0.95)
    
    Q = qunitroot(X, 100, trend = "c",   statistic = "t")
    P = punitroot(Q, 100, trend = "c",   statistic = "t")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)
    
    Q = qunitroot(X, 100, trend = "nc",  statistic = "t")
    P = punitroot(Q, 100, trend = "nc",  statistic = "t")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)
    
    Q = qunitroot(X, 100, trend = "ct",  statistic = "t")
    P = punitroot(Q, 100, trend = "ct",  statistic = "t")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)
    
    Q = qunitroot(X, 100, trend = "ctt", statistic = "t")
    P = punitroot(Q, 100, trend = "ctt", statistic = "t")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)
    
    Q = qunitroot(X, 100, trend = "c",   statistic = "n")
    P = punitroot(Q, 100, trend = "c",   statistic = "n")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)
    
    Q = qunitroot(X, 100, trend = "nc",  statistic = "n")
    P = punitroot(Q, 100, trend = "nc",  statistic = "n")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)
    
    Q = qunitroot(X, 100, trend = "ct",  statistic = "n")
    P = punitroot(Q, 100, trend = "ct",  statistic = "n")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol) 
    
    Q = qunitroot(X, 100, trend = "ctt", statistic = "n")
    P = punitroot(Q, 100, trend = "ctt", statistic = "n")
    print(cbind(Q, P))
    checkEqualsNumeric(target = X, current = P, tolerance = tol)

    # Return Value:
    return()    
}


################################################################################

