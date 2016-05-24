
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
# FUNCTION:                     BINARY OPTIONS:
#  GapOption                     Gap Option
#  CashOrNothingOption           Cash Or Nothing Option
#  TwoAssetCashOrNothingOption   Two Asset Cash-Or Nothing Option
#  AssetOrNothingOption          Asset Or Nothing Option
#  SuperShareOption              Super Share Option
#  BinaryBarrierOption           Binary Barrier Option
################################################################################


test.GapOption = 
function()
{
    # Examples from Chapter 2.11 in E.G. Haug's Option Guide (1997)
    
    # Gap Option [2.11.1]:
    GapOption(TypeFlag = "c", S = 50, X1 = 50, X2 = 57, Time = 0.5, 
        r = 0.09, b = 0.09, sigma = 0.20)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.CashOrNothingOption = 
function()
{    
    # Cash Or Nothing Option [2.11.2]:
    CashOrNothingOption(TypeFlag = "p", S = 100, X = 80, K = 10, 
        Time = 9/12, r = 0.06, b = 0, sigma = 0.35) 
    
    # Two Asset Cash Or Nothing Option [2.11.3]:
    # Type 1 - call:
    TwoAssetCashOrNothingOption(TypeFlag = "c", S1 = 100, S2 = 100, 
        X1 = 110, X2 = 90, K = 10, Time = 0.5, r = 0.10, b1 = 0.05, 
        b2 = 0.06, sigma1 = 0.20, sigma2 = 0.25, rho = 0.5)
    # Type 2 - put:
    TwoAssetCashOrNothingOption(TypeFlag = "p", S1 = 100, S2 = 100, 
        X1 = 110, X2 = 90, K = 10, Time = 0.5, r = 0.10, b1 = 0.05, 
        b2 = 0.06, sigma1 = 0.20, sigma2 = 0.25, rho = -0.5)
    # Type 3 - down-up:
    TwoAssetCashOrNothingOption(TypeFlag = "ud", S1 = 100, S2 = 100, 
        X1 = 110, X2 = 90, K = 10, Time = 1, r = 0.10, b1 = 0.05, 
        b2 = 0.06, sigma1 = 0.20, sigma2 = 0.25, rho = 0)
    # Type 4 - up-down:
    TwoAssetCashOrNothingOption(TypeFlag = "du", S1 = 100, S2 = 100, 
        X1 = 110, X2 = 90, K = 10, Time = 1, r = 0.10, b1 = 0.05, 
        b2 = 0.06, sigma1 = 0.20, sigma2 = 0.25, rho = 0)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.AssetOrNothingOption = 
function()
{    
    # Asset Or Nothing Option [2.11.4]: 
    AssetOrNothingOption(TypeFlag = "p", S = 70, X = 65, Time = 0.5, 
        r = 0.07, b = 0.07 - 0.05, sigma = 0.27)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.SuperShareOption = 
function()
{    
    # Super Share Option [2.11.5]:  
    SuperShareOption(S = 100, XL = 90, XH = 110, Time = 0.25, r = 0.10, 
        b = 0, sigma = 0.20)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.BinaryBarrierOption = 
function()
{    
    # Binary Barrier Option [2.11.6]: 
    BinaryBarrierOption(TypeFlag = "6", S = 95, X=102, H = 100, 
        K = 15, Time = 0.5, r = 0.1, b = 0.1, sigma = 0.20)
    BinaryBarrierOption(TypeFlag = "12", S = 95, X = 98, H = 100, 
        K = 15, Time = 0.5, r = 0.1, b = 0.1, sigma = 0.20)
        
    # Return Value:
    return()    
}


################################################################################

