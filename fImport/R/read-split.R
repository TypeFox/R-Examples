
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

# Copyrights (C) for this R-port:
#   1999 - 2012 Diethelm Wuertz, Zurich, <wuertz@itp.phys.ethz.ch>
#   2009 - 2012 Rmetrics Association, Zurich, www.rmetrics.org


################################################################################
# FUNCTION:               DESCRIPTION:
#  dataSplit               Splits data matrix from a downloaded file
#  charvecSplit            Splits charvec vector from a downloaded file
################################################################################


dataSplit <- 
function (x, split=" ", col=-1) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description
    #   Splits data matrix from a downloaded file
    
    # FUNCTION
    
    # Split Function:
    FUN = function(x, split, col) {
        unlist(strsplit(x, split))[col]
    }
    
    # Split Data:
    data = unlist(lapply(x, FUN, split = split, col = col))
    data = matrix(as.numeric(data), byrow = TRUE, nrow = length(x))
    
    # Return Value:
    data
}


# ------------------------------------------------------------------------------


charvecSplit <-
function (x, split=" ", col=1, format="%F") 
{
    # A function implemented by Diethelm Wuertz
    
    # Description
    #   Splits charvec vector from a downloaded file
    
    # FUNCTION
    
    # Split Function:
    FUN = function(x, split, col) {
        unlist(strsplit(x, split))[col]
    }
    
    # Split Date Character Vector:
    charvec = unlist(lapply(x, FUN, split = split, col = col))
    charvec = format(strptime(charvec, format = format))
    
    # Return Value:
    charvec
}


################################################################################

