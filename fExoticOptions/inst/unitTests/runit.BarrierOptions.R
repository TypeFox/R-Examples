
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
# FUNCTION:                     BARRIER OPTIONS:
#  StandardBarrierOption         Standard Barrier Option
#  DoubleBarrierOption           Double Barrier Option
#  PTSingleAssetBarrierOption    Partial Time Barrier Option
#  TwoAssetBarrierOption         Two Asset Barrier
#  PTTwoAssetBarrierOption       Partial Time TwoAsset Barrier Option
#  LookBarrierOption             Look Barrier Option
#  DiscreteBarrierOption         Discrete Adjusted Barrier Option
#  SoftBarrierOption             Soft Barrier Option
################################################################################


test.StandardBarrierOption = 
function()
{
    # Examples from Chapter 2.10 in E.G. Haug's Option Guide (1997)
    
    # Standard Barrier Option [2.10.1]:
    # down-and-out Barrier Call
    StandardBarrierOption(TypeFlag = "cdo", S = 100, X = 90, 
        H = 95, K = 3, Time = 0.5, r = 0.08, b = 0.04, sigma = 0.25)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.DoubleBarrierOption = 
function()
{    
    # Double Barrier Option [2.10.2]:
    DoubleBarrierOption(TypeFlag = "co", S = 100, X = 100, L = 50, 
        U = 150, Time = 0.25, r = 0.10, b = 0.10, sigma = 0.15, 
        delta1 = -0.1, delta2 = 0.1)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.PTSingleAssetBarrierOption = 
function()
{    
    # Partial Time Single-Asset Barrier Option [2.10.3]:
    PTSingleAssetBarrierOption(TypeFlag = "coB1", S = 95, X = 110, 
        H = 100, time1 = 0.5, Time2 = 1, r = 0.20, b = 0.20, 
        sigma = 0.25)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.TwoAssetBarrierOption = 
function()
{    
    # Two Asset Barrier Option [2.10.4]:
    TwoAssetBarrierOption(TypeFlag = "puo", S1 = 100, S2 = 100, 
        X = 110, H = 105, Time = 0.5, r = 0.08, b1 = 0.08, b2 = 0.08, 
        sigma1 = 0.2, sigma2 = 0.2, rho = -0.5)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.PTTwoAssetBarrierOption = 
function()
{    
    # PT Two Asset Barrier Option [2.10.5]:
    PTTwoAssetBarrierOption(TypeFlag = "pdo", S1 = 100, S2 = 100, 
        X = 100, H = 85, time1 = 0.5, Time2 = 1, r = 0.1, b1 = 0.1, 
        b2 = 0.1, sigma1 = 0.25, sigma2 = 0.30, rho = -0.5)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.LookBarrierOption = 
function()
{    
    # Look Barrier Option [2.10.6]:
    LookBarrierOption(TypeFlag = "cuo", S = 100, X = 100, H = 130, 
        time1 = 0.25, Time2 = 1, r = 0.1, b = 0.1, sigma = 0.15)
    LookBarrierOption(TypeFlag = "cuo", S = 100, X = 100, H = 110, 
        time1 = 1, Time2 = 1, r = 0.1, b = 0.1, sigma = 0.30)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.DiscreteBarrierOption = 
function()
{    
    # Discrete Barrier Option [2.10.7]:  
    DiscreteBarrierOption(S = 100, H = 105, sigma = 0.25, dt = 0.1)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.SoftBarrierOption = 
function()
{    
    # Soft Barrier Option [2.10.8]:
    SoftBarrierOption(TypeFlag = "cdo", S = 100, X = 100, L = 70, 
        U = 95, Time = 0.5, r = 0.1, b = 0.05, sigma = 0.20)

    # Return Value:
    return()    
}


################################################################################

