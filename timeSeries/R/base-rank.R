#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  ../../COPYING


################################################################################
# FUNCTION:                 COLUMN STATISTICS IN FUTILITIES:
#  rank,timeSeries           Returns sample ranks of a 'timeSeries' object
################################################################################


setMethod("rank", "timeSeries",
    function(x,  na.last = TRUE,
        ties.method = eval(formals(rank)$ties.method))
    {
        # Description:
        #   Returns the sample ranks of the values in a 'timeSeries' 
        
        # Arguments:
        #   x - an object of class 'timeSeries'
        #   ties.method - 
        #       "average", replaces them by their mean, 
        #       "first" method results in a permutation with increasing 
        #           values at each index set of ties. 
        #       "random" method puts these in random order whereas the 
        #           default, 
        #       "max" and "min" replaces them by their maximum and minimum 
        #           respectively, the latter being the typical sports ranking. 

        # Note:
        #   Ties (i.e., equal values) and missing values can be handled 
        #   in several ways. 

        # FUNCION:
        
        # Return Value:
        apply(x, 2, rank, na.last = na.last, ties.method = ties.method)
    }
)


################################################################################

