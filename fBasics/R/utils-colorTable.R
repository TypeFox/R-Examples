
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
# You should have received A copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA  02111-1307  USA


################################################################################
# FUNCTION:                 DESCRIPTION:
#  colorTable                Shows a table of plot color codes
################################################################################


colorTable <- 
function(cex = 0.7) 
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Displays a table of plot colors.
    
    # Author:
    #   Unknown, piece of code found on the internet.
    
    # Example:
    #   colorTable()
  
    # FUNCTION:
    
    # Plot:
    plot(0, 0, xlim = c(-1, 10), ylim = c(0, 10), type = 'n', axes = FALSE, 
        xlab = '', ylab = '', cex = cex, main = "Table of Color Codes")
    ## DW
    ## j = -1
    j = -2 ## here I have shifted the plot to the left ... DW
    for(i in 0:99) {
        if(i %% 10 == 0) {j = j+1; k = 10}
        k = k-1
        points(j, k, pch = 15, col = i, cex = 2)
        text(j+0.45, k, i, cex = cex)}
    
    # Return Value:
    invisible()
}


################################################################################

