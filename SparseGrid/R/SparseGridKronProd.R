# Copyright (C) 2011 Jelmer Ypma. All Rights Reserved.
# This code is published under the GPL.
#
# File:   SparseGridKronProd.R
# Author: Jelmer Ypma
# Date:   31 March 2011
#
# This code is based on code from www.sparse-grids.de
# with permission from the authors.
#
# Description: function for generating tensor product quadrature rule.
#			   for internal use.
#
# Input: 
#     n1D = cell array of 1D nodes (nodes are vectors, otherwise length doesn't work)
#     w1D = cell array of 1D weights (weights are vectors, oterwise length doesn't work)
#
# Output: 
#	list with two elements:
#     nodes   = matrix of tensor product nodes 
#     weights = vector of tensor product weights 
#

SparseGridKronProd <- function( n1D, w1D ) {
    nodes   <- matrix( n1D[[ 1 ]], nrow=length( n1D[[ 1 ]] ), ncol=1 )
    weights <- w1D[[ 1 ]]

    # check if length > 1, because otherwise 2:length would be c(2,1)
    if ( length( n1D ) > 1 ) {
        for ( j in 2:length(n1D) ) {
            newnodes    <- n1D[[ j ]]
            nodes       <- cbind( kronecker( nodes, rep(1, length(newnodes))), 
                                  kronecker( rep(1, nrow(nodes)), newnodes) )
            weights     <- kronecker( weights, w1D[[ j ]] )
        }
    }
    
    return( list( "nodes"     = nodes,
                  "weights"   = weights ) )
}
