
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
# FUNCTION:                           LOOKBACK OPTIONS:
#  FloatingStrikeLookbackOption        Floating Strike Lookback Option
#  FixedStrikeLookbackOption           Fixed Strike Lookback Option
#  PTFloatingStrikeLookbackOption      Partial Floating Strike LB Option
#  PTFixedStrikeLookbackOption         Partial Fixed Strike LB Option  
#  ExtremeSpreadOption                 Extreme Spread Option
################################################################################


test.FloatingStrikeLookbackOption = 
function()
{
    # Examples from Chapter 2.9 in E.G. Haug's Option Guide (1997)

    # Floating Strike Lookback Option [2.9.1]:
    FloatingStrikeLookbackOption(TypeFlag = "c", S = 120, 
        SMinOrMax = 100, Time = 0.5, r = 0.10, b = 0.10-0.06, 
        sigma = 0.30)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.FixedStrikeLookbackOption = 
function()
{    
    # Fixed Strike Lookback Option [2.9.2]:
    FixedStrikeLookbackOption(TypeFlag = "c", S = 100, 
        SMinOrMax = 100, X = 105, Time = 0.5, r = 0.10, b = 0.10, 
        sigma = 0.30)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.PTFloatingStrikeLookbackOption = 
function()
{       
    # Partial Time Floating Strike Lookback Option [2.9.3]:
    PTFloatingStrikeLookbackOption(TypeFlag = "p", S = 90, 
        SMinOrMax = 90, time1 = 0.5, Time2 = 1, r = 0.06, b = 0.06, 
        sigma = 0.20, lambda  = 1)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.PTFixedStrikeLookbackOption = 
function()
{       
    # Partial Time Fixed Strike Lookback Option [2.9.4]:
    PTFixedStrikeLookbackOption(TypeFlag = "c", S = 100, X = 90, 
        time1 = 0.5, Time2 = 1, r = 0.06, b = 0.06, sigma = 0.20)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.ExtremeSpreadOption = 
function()
{         
    # Extreme Spread Option [2.9.5]:
    ExtremeSpreadOption(TypeFlag = "c", S = 100, SMin = NA, 
        SMax = 110, time1 = 0.5, Time2 = 1, r = 0.1, b = 0.1, 
        sigma = 0.30)
    ExtremeSpreadOption(TypeFlag = "cr", S = 100, SMin = 90, 
        SMax = NA, time1 = 0.5, Time2 = 1, r = 0.1, b = 0.1, 
        sigma = 0.30)

    # Return Value:
    return()    
}


################################################################################

