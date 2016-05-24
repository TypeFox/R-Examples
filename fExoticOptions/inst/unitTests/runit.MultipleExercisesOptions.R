
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
# FUNCTION:                       MULTIPLE EXERCISES OPTIONS:           
#  ExecutiveStockOption            Executive Stock Option
#  ForwardStartOption              Forward Start Option
#  RatchetOption                   Ratchet [Compound] Option
#  TimeSwitchOption                Time Switch Option
#  SimpleChooserOption             Simple Chooser Option
#  ComplexChooserOption            Complex Chooser Option
#  OptionOnOption                  Options On Options
#  HolderExtendibleOption          Holder Extendible Option
#  WriterExtendibleOption          Writer Extendible Option
################################################################################


test.ExecutiveStockOption = 
function()
{
    # Examples from Chapter 2.1 - 2.7 in E.G. Haug's Option Guide (1997)
    
    # ExecutiveStockOption [2.1]:
    ExecutiveStockOption(TypeFlag = "c", S = 60, X = 64, Time = 2, 
        r = 0.07, b = 0.07-0.03, sigma = 0.38, lambda = 0.15) 
      
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.ForwardStartOption = 
function()
{   
    # ForwardStartOption [2.2]:
    ForwardStartOption(TypeFlag = "c", S = 60, alpha = 1.1, 
        time1 = 1, Time2 = 1/4, r = 0.08, b = 0.08-0.04, sigma = 0.30)
      
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.RatchetOption = 
function()
{     
    # Ratchet Option [2.3]:
    RatchetOption(TypeFlag = "c", S = 60, alpha = 1.1, time1 = c(1.00, 0.75), 
        Time2 = c(0.75, 0.50), r = 0.08, b = 0.04, sigma = 0.30)
      
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.TimeSwitchOption = 
function()
{     
    # Time Switch Option [2.4]:
    TimeSwitchOption(TypeFlag = "c", S = 100, X = 110, Time = 1, 
        r = 0.06, b = 0.06, sigma = 0.26, A = 5, m = 0, dt = 1/365)
      
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.SimpleChooserOption = 
function()
{    
    # Simple Chooser Option [2.5.1]:
    SimpleChooserOption(S = 50, X = 50, time1 = 1/4, Time2 = 1/2, 
        r = 0.08, b = 0.08, sigma = 0.25)  
      
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.ComplexChooserOption = 
function()
{        
    # Complex Chooser Option [2.5.2]:
    ComplexChooserOption(S = 50, Xc = 55, Xp = 48, Time = 0.25, 
        Timec = 0.50, Timep = 0.5833, r = 0.10, b = 0.1-0.05, 
        sigma = 0.35, doprint = TRUE)
      
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.OptionOnOption = 
function()
{     
    # Option On Option [2.6]:
    OptionOnOption(TypeFlag = "pc", S = 500, X1 = 520, X2 = 50, 
        time1 = 1/2, Time2 = 1/4, r = 0.08, b = 0.08-0.03, sigma = 0.35)
      
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.HolderExtendibleOption = 
function()
{    
    # Holder Extendible Option [2.7.1]:
    HolderExtendibleOption(TypeFlag = "c", S = 100, X1 = 100, 
        X2 = 105, time1 = 0.50, Time2 = 0.75, r = 0.08, b = 0.08, 
        sigma = 0.25, A = 1)
      
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.WriterExtendibleOption = 
function()
{
    # Writer Extendible Option [2.7.2]:
    WriterExtendibleOption(TypeFlag = "c", S = 80, X1 = 90, X2 = 82,
        time1 = 0.50, Time2 = 0.75, r = 0.10, b = 0.10, sigma = 0.30)
    
    # Return Value:
    return()    
}


################################################################################

