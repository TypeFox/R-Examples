
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA 02111-1307 USA


################################################################################
# FUNCTION:                     DESCRIPTION:
#  'fPORTFOLIOCON'               S4 Portfolio Constraints Class
################################################################################


setClass("fPFOLIOCON", 
    representation(
    
        # Function Implemented by Diethelm Wuertz
    
        # Constraints expressed by strings
        #   LongOnly, Short, Full, ...
        stringConstraints = "character",
        
        # BoxConstraints:
        minWConstraints = "numeric",
        maxWConstraints = "numeric",
        
        # Group Constraints:
        eqsumWConstraints = "matrix",
        minsumWConstraints = "matrix",
        maxsumWConstraints = "matrix",
        
        # Covariance Risk Budget Constraints:
        minBConstraints = "numeric", 
        maxBConstraints = "numeric",
        
        # Nonlinear Constraints:
        listFConstraints = "list",
        minFConstraints = "numeric", 
        maxFConstraints = "numeric",
        
        # Buyin Constraints:
        minBuyinConstraints = "numeric", 
        maxBuyinConstraints = "numeric",
        
        # Cardinality Constraints:
        nCardConstraints = "integer",
        minCardConstraints = "numeric", 
        maxCardConstraints = "numeric")
        
) 
        

################################################################################

