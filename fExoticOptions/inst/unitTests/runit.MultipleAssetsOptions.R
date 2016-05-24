
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
# FUNCTION:                       MULTI ASSET OPTION:
#  TwoAssetCorrelationOption       Two Asset Correlation Option
#   [ExchangeOneForAnotherOption]    [Exchange One For Another Option]  
#  EuropeanExchangeOption          European Exchange Optionn
#  AmericanExchangeOption          American Exchange Option
#  ExchangeOnExchangeOption        Exchange Exchange Option
#  TwoRiskyAssetsOption            Option On The MinMax
#  SpreadApproxOption              Spread Approximated Option              
################################################################################


test.TwoAssetCorrelationOption = 
function()
{
    # Examples from Chapter 2.8 in E.G. Haug's Option Guide (1997)
    
    # Two Asset Correlation Options [2.8.1]:
    TwoAssetCorrelationOption(TypeFlag = "c", S1 = 52, S2 = 65, 
        X1 = 50, X2 = 70, Time = 0.5, r = 0.10, b1 = 0.10, b2 = 0.10, 
        sigma1 = 0.2, sigma2 = 0.3, rho = 0.75) 

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.EuropeanExchangeOption = 
function()
{    
    # European Exchange Options [2.8.2]: 
    EuropeanExchangeOption(S1 = 22, S2 = 0.20, Q1 = 1, Q2 = 1, 
        Time = 0.1, r = 0.1, b1 = 0.04, b2 = 0.06, sigma1 = 0.2, 
        sigma2 = 0.25, rho = -0.5)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.AmericanExchangeOption = 
function()
{    
    # American Exchange Options [2.8.2]:
    AmericanExchangeOption(S1 = 22, S2 = 0.20, Q1 = 1, Q2 = 1, 
        Time = 0.1, r = 0.1, b1 = 0.04, b2 = 0.06, sigma1 = 0.2, 
        sigma2 = 0.25, rho = -0.5)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.ExchangeOnExchangeOption = 
function()
{    
    # Exchange Options On Exchange Options [2.8.3]:
    for (flag in 1:4) print(
    ExchangeOnExchangeOption(TypeFlag = as.character(flag), 
        S1 = 105, S2 = 100, Q = 0.1, time1 = 0.75, Time2 = 1.0, r = 0.1, 
        b1 = 0.10, b2 = 0.10, sigma1 = 0.20, sigma2 = 0.25, rho = -0.5))

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.TwoRiskyAssetsOption = 
function()
{    
    # Two Risky Assets Options [2.8.4]:
    TwoRiskyAssetsOption(TypeFlag = "cmax", S1 = 100, S2 = 105, 
        X = 98, Time = 0.5, r = 0.05, b1 = -0.01, b2 = -0.04, 
        sigma1 = 0.11, sigma2 = 0.16, rho = 0.63)
    TwoRiskyAssetsOption(TypeFlag = "pmax", S1 = 100, S2 = 105, 
        X = 98, Time = 0.5, r = 0.05, b1 = -0.01, b2 = -0.04, 
        sigma1 = 0.11, sigma2 = 0.16, rho = 0.63)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.SpreadApproxOption = 
function()
{    
    # Spread-Option Approximation [2.8.5]:
    SpreadApproxOption(TypeFlag = "c", S1 = 28, S2 = 20, X = 7, 
        Time = 0.25, r = 0.05, sigma1 = 0.29, sigma2 = 0.36, rho = 0.42)
         
    # Return Value:
    return()    
}


################################################################################

