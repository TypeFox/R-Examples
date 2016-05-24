
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
#  RollGeskeWhaleyOption      Roll-Geske-Whaley Calls on Dividend Paying Stocks
#  BAWAmericanApproxOption    Barone-Adesi and Whaley Approximation
#  BSAmericanApproxOption     Bjerksund and Stensland Approximation
################################################################################


test.RollGeskeWhaleyOption =
function()
{
    # RollGeskeWhaleyOption      
    #   Roll-Geske-Whaley Calls on Dividend Paying Stocks
    
    # Arguments:
    #   RollGeskeWhaleyOption(S, X, time1, Time2, r, D, sigma, 
    #   title = NULL, description = NULL)
    
    # Roll-Geske-Whaley American Calls on Dividend Paying 
    # Stocks [Haug 1.4.1]
    RollGeskeWhaleyOption(S = 80, X = 82, time1 = 1/4, 
        Time2 = 1/3, r = 0.06, D = 4, sigma = 0.30)
   
    # Return Value:
    return()    
}   


# ------------------------------------------------------------------------------


test.BAWAmericanApproxOption =
function()
{
    # BAWAmericanApproxOption    
    #   Barone-Adesi and Whaley Approximation
      
    # Arguments:
    #   BAWAmericanApproxOption(TypeFlag = c("c", "p"), S, X, Time, r, b, 
    #   sigma, title = NULL, description = NULL)
    
    # Barone-Adesi and Whaley Approximation for American 
    # Options [Haug 1.4.2] vs. Black76 Option on Futures:
    BAWAmericanApproxOption(TypeFlag = "p", S = 100, 
        X = 100, Time = 0.5, r = 0.10, b = 0, sigma = 0.25)
    Black76Option(TypeFlag = "c", FT = 100, X = 100, 
        Time = 0.5, r = 0.10, sigma = 0.25)  
  
    # Return Value:
    return()    
}    


# ------------------------------------------------------------------------------


test.BSAmericanApproxOption =
function()
{
    # BSAmericanApproxOption     
    #   Bjerksund and Stensland Approximation  
     
    # Arguments:
    #   BSAmericanApproxOption(TypeFlag = c("c", "p"), S, X, Time, r, b, 
    #   sigma, title = NULL, description = NULL) 
   
    # Bjerksund and Stensland Approximation for American Options:
    BSAmericanApproxOption(TypeFlag = "c", S = 42, X = 40, 
        Time = 0.75, r = 0.04, b = 0.04-0.08, sigma = 0.35)
     
    # Return Value:
    return()    
}     


################################################################################

