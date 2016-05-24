
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
# FUNCTION:                     CURRENCY TRANSLATED OPTIONS:
#  FEInDomesticFXOption          FX In Domestic Currency
#  QuantoOption                  Quanto Option
#  EquityLinkedFXOption          EquityLinked FX Option
#  TakeoverFXOption              Takeover FX Option
################################################################################


test.FEInDomesticFXOption = 
function()
{
    # Examples from Chapter 2.13 in E.G. Haug's Option Guide (1997)
    
    # Foreign Equity Options Struck in Domestic Currency [2.13.1]:
    FEInDomesticFXOption(TypeFlag = "c", S = 100, E = 1.5, 
        X = 160, Time = 0.5, r = 0.08, q = 0.05, sigmaS = 0.20, 
        sigmaE = 0.12, rho = 0.45)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.QuantoOption = 
function()
{    
    # Fixed Exchange-Rate Foreign-Equity Option [2.13.2]: 
    QuantoOption(TypeFlag = "c", S = 100, Ep = 1.5, X = 105, 
        Time = 0.5, r = 0.08, rf = 0.05, q = 0.04, sigmaS= 0.2, 
        sigmaE = 0.10, rho = 0.30) 

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.EquityLinkedFXOption = 
function()
{    
    # Equity Linked Foreign Exchange Option [2.13.3]:
    EquityLinkedFXOption(TypeFlag = "p", E = 1.5, S = 100, 
        X = 1.52, Time = 0.25, r = 0.08, rf = 0.05, q = 0.04, 
        sigmaS = 0.20, sigmaE = 0.12, rho = -0.40)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.TakeoverFXOption = 
function()
{    
    # Takeover Foreign-Exchange Option [2.13.4]:
    TakeoverFXOption(V = 100, B = 100, E = 1.5, X = 1.55, Time = 1, 
        r = 0.08, rf = 0.06, sigmaV = 0.20, sigmaE = 0.25, rho = 0.1)

    # Return Value:
    return()    
}


################################################################################

