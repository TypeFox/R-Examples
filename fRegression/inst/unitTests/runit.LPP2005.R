
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


###########################\#####################################################


test.Fit <- 
    function()
{   
    # Simulate Artificial LM:
    x <- as.timeSeries(data(LPP2005REC))
    
    # Fit Parameters:
    lmfit <- regFit(LPP40 ~ 0 + SPI + SBI + SII + LMI + MPI + ALT, 
        data = x, use = "lm") 
    
    # 
    print(lmfit)
 
    # Return Value:
    return()
}


###############################################################################


