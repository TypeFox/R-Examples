
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


################################################################################
# FUNCTION:           DESCRIPTION:
#  squareBinning       Square binning of irregularly distributed data sets
#  plot                S3 Method for plotting square binned data sets
# FUNCTION:           DESCRIPTION:
#  hexBinning          Hexagonal binning of irregularly distributed data sets
#  plot                S3 Method for plotting hexagonal binned data sets
################################################################################


test.squareBinning = 
function()
{
    #  squareBinning    Square binning of irregularly distributed data sets
    #  plot             S3 Method for plotting square binned data sets

    # Generate Grid Data:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    sB = squareBinning(x = rnorm(1000), y = rnorm(1000))
    
    # Plot:
    par(mfrow = c(1, 1))
    plot(sB)
    title(main = "Square Binning")
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.hexBinning = 
function()
{
    #  hexBinning       Hexagonal binning of irregularly distributed data sets
    #  plot             S3 Method for plotting hexagonal binned data sets
    
    # Generate Grid Data:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    hB = hexBinning(x = rnorm(1000), y = rnorm(1000))
    
    # Plot:
    par(mfrow = c(1, 1))
    plot(hB)
    title(main = "Hexagonal Binning")
      
    # Return Value:
    return()    
}


################################################################################

