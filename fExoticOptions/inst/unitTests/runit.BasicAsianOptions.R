
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
# FUNCTION:                         ASIAN OPTIONS:
#  GeometricAverageRateOption        Geometric Average Rate Option
#  TurnbullWakemanAsianApproxOption  Turnbull-Wakeman Approximated Asian Option
#  LevyAsianApproxOption             Levy Approximated Asian Option
################################################################################


test.GeometricAverageRateOption = 
function()
{
    # Examples from Chapter 2.12 in E.G. Haug's Option Guide (1997)
    
    # Geometric Average Rate Option:
    GeometricAverageRateOption(TypeFlag = "p", S = 80, X = 85, 
        Time = 0.25, r = 0.05, b = 0.08, sigma = 0.20)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.TurnbullWakemanAsianApproxOption = 
function()
{    
    # Turnbull Wakeman Approximation:
    TurnbullWakemanAsianApproxOption(TypeFlag = "p", S = 90, SA = 88, 
        X = 95, Time = 0.50, time = 0.25, tau = 0.0, r = 0.07, 
        b = 0.02, sigma = 0.25)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.LevyAsianApproxOption = 
function()
{    
    # Levy Asian Approximation:   
    LevyAsianApproxOption(TypeFlag = "c", S = 100, SA = 100, X = 105, 
        Time = 0.75, time = 0.50, r = 0.10, b = 0.05, sigma = 0.15)
 
    # Return Value:
    return()    
}


################################################################################

