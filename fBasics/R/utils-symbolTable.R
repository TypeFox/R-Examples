
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
#  symbolTable               Shows a table of plot symbols from a given font
################################################################################


symbolTable <- 
function(font = par('font'), cex = 0.7) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Shows a table of plot characters from a given font
    
    # Example:
    #   symbolTable()
    
    # Author:
    #   Unknown, piece of code found on the internet.

    # FUNCTION:
    
    # Table:
    plot(0, 0, xlim = c(-1, 11), ylim = c(0, 26), type = 'n', 
        axes = FALSE, xlab = '', ylab = '', 
        main = "Table of Plot Characters")
    j = -1
    for(i in 0:255) {
        if(i %% 25 == 0) {j = j+1; k = 26}
        k = k-1
        points(j, k, pch = i, font = font, cex = cex, col = 2)
        text(j + 0.50, k, i, cex = cex) 
    }
    
    # Return Value:
    invisible(font)
}


################################################################################

