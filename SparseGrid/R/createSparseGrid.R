# Copyright (C) 2011 Jelmer Ypma. All Rights Reserved.
# This code is published under the GPL.
#
# File:   createSparseGrid.R
# Author: Jelmer Ypma
# Date:   31 March 2011
#
# This code is based on code from www.sparse-grids.de
# with permission from the authors.
#
# Description: function to calculate nodes and weights for numerical integration 
#              on sparse grids (Smolyak).
#
# Input: 
#     type       : String for type of 1D integration rule:
#                  "KPU": Nested rule for unweighted integral over [0,1]
#                  "KPN": Nested rule for integral with Gaussian weight
#                  "GQU": Gaussian quadrature for unweighted integral over [0,1] (Gauss-Legendre)
#                  "GQN": Gaussian quadrature for integral with Gaussian weight (Gauss-Hermite)
#                  func:  any function name. Function must accept level l and
#                         return nodes n and weights w for univariate quadrature rule with
#                         polynomial exactness 2l-1 as [n w] = feval(func,level)
#     dimension  : dimension of the integration problem
#     k          : Accuracy level. The rule will be exact for polynomial up to total order 2k-1
#     sym        : (optional) only used for own 1D quadrature rule (type not "KPU",...). If
#                  sym is supplied and not=0, the code will run faster but will
#                  produce incorrect results if 1D quadrature rule is asymmetric
#
# Output: 
#     list with two elements:
#		nodes    	= matrix of nodes with dim columns 
#     	weights    = row vector of corresponding weights
#
#
# 14/05/2012: Fixed a bug related to checking of arguments when 'type' is a user-supplied function.
#
# TODO: check for sum(weights)==1 and length(weights)==nrow(nodes)

createSparseGrid <- function(type, dimension, k, sym = FALSE) {

    # assertions    
	if (! as.integer(dimension) == dimension) {
		stop( paste("createSparseGrid expects dimension to be integer, user supplied dimension = ", dimension, sep='') )
	}
	if (! as.integer(k) == k) {
		stop( paste("createSparseGrid expects k to be integer, user supplied k = ", k, sep='') )
	}
	if (! is.logical(sym) | is.na(sym)) {
		stop( paste("createSparseGrid expects sym to be logical, user supplied sym = ", sym, sep='') )
	}
    
    # check if a function or a string was supplied
    # builtinfct is a boolean defining whether 'type' is a built-in function (TRUE) or a user-supplied function (FALSE)
    builtinfct <- FALSE
	if ( is.function(type) ) {
		if (length( formals( type ) ) != 1 ) {
			stop( "User supplied function (argument type) to createSparseGrid needs to have one argument." )
		}
		tmp.grid <- type( dimension )
		if (! is.list(tmp.grid) ) {
			stop( "User supplied function (argument type) to createSparseGrid needs to return a list." )
		}
		if (! is.list(tmp.grid) ) {
			stop( "User supplied function (argument type) to createSparseGrid needs to return a list." )
		}
		if (! "nodes" %in% names(tmp.grid) ) {
			stop( "Element 'nodes' not found in list returned by user supplied function (argument type) to createSparseGrid." )	
		}
		if (! "weights" %in% names(tmp.grid) ) {
			stop( "Element 'weights' not found in list returned by user supplied function (argument type) to createSparseGrid." )	
		}
	}
    else if ( is.character(type) ) {
        # input string is from the list of built-in functions?
        builtinfct <- type %in% c('GQU', 'GQN', 'KPU', 'KPN' )
        if ( builtinfct ) { 
            sym = TRUE
        } else {
            stop( paste("createSparseGrid expects type to be a string with values 'GQU', 'GQN', 'KPU', 'KPN', or a function, user supplied type = ", type, sep='') )
        }
    }
    else {
        stop( paste("createSparseGrid expects type to be a string with values 'GQU', 'GQN', 'KPU', 'KPN', or a function, user supplied type = ", type, sep='') )
    }
	
	if ( !builtinfct & sym ) {
		warning('Using sym=TRUE for non-symmetric rules produces incorrect results. Rules that are symmetric in theory can be asymmetric in practice because of numerical errors. See demo_userdefined.R for an example where this can go wrong.')
	}
    
    # get 1D nodes & weights
    tryCatch( {
        # create empty lists of length k
        n1D <- vector( mode="list", length=k )
        w1D <- vector( mode="list", length=k )
        R1D <- rep( 0, k )      # vector with zeros
        for ( level in 1:k ) {
            if ( builtinfct ) {
                # 'type' is a string
                res         <- eval(call(type,level))
            } else {
                # 'type' is a function
                res         <- type( level )
            }
            nodes       <- res$nodes 
            weights     <- res$weights
            
            # if user supplied symmetric 1D rule: keep only positive orthant
            if ( !builtinfct & sym ) {
                num.new         <- length( nodes )
                nodes.sorted    <- sort( nodes, index.return=TRUE )    # sortvec
                nodes           <- nodes.sorted$x
                weights         <- weights[ nodes.sorted$ix ]
                nodes           <- nodes[ (floor(num.new/2)+1):num.new ]
                weights         <- weights[ (floor(num.new/2)+1):num.new ]
            }
            R1D[  level ]   <- length( weights )
            n1D[[ level ]]  <- nodes
            w1D[[ level ]]  <- weights
        }
    }, error = function(e) { cat("Error evaluating the 1D rule\n") } )
    
    
    # initialization
    minq    <- max( 0, k - dimension )
    maxq    <- k-1
    nodes   <- matrix(0, nrow=0, ncol=dimension)
    weights <- numeric(0)
    
    # outer loop over q
    for ( q.cnt in minq:maxq ) {
        r   <- length(weights)
        bq  <- (-1)^(maxq - q.cnt) * choose( dimension - 1, dimension + q.cnt - k )
        
        # matrix of all rowvectors in N^D_{q.cnt}
        indices.mat = SparseGridGetSeq( dimension, dimension + q.cnt )
        
        # preallocate new rows for nodes & weights
        Rq      <- sapply( 1:nrow(indices.mat), function(row.cnt) { prod(R1D[indices.mat[row.cnt,]]) } )
        Rq.sum  <- sum( Rq )
        nodes   <- rbind( nodes, matrix(0, nrow=Rq.sum, ncol=dimension) )
        weights <- c( weights, rep(0, Rq.sum) )
        
        # inner loop collecting product rules
        for ( j in 1:nrow(indices.mat) ) {
            midx        <- indices.mat[ j, ]
            res         <- SparseGridKronProd( n1D[midx],w1D[midx] )
            nodes[ (r+1):(r+Rq[j]), ]   <- res$nodes
            weights[ (r+1):(r+Rq[j]) ]  <- bq * res$weights        # .*
            r = r + Rq[j]
        }
        
        # collect equal nodes: first sort
        nodes.sorted    <- sortrows( nodes, index.return=TRUE )    # sortvec
        nodes           <- nodes.sorted$x
        weights         <- weights[ nodes.sorted$ix ]
        keep            <- 1 
        lastkeep        <- 1
        
        # then make a list of rows to keep and sum weights of equal nodes
        if ( nrow(nodes) > 1 ) {
            for ( j in 2:nrow( nodes ) ) {
                if ( all( nodes[ j, ] == nodes[ j-1, ] ) ) {
                    weights[lastkeep] <- weights[ lastkeep ] + weights[ j ]
                } 
                else {
                  lastkeep  <- j
                  keep      <- c( keep, j )
                }
            }
        }
        
        # drop equal rows
        # (add matrix statement to keep nodes as a matrix, if it contains only 1 column)
        nodes   <- matrix( nodes[ keep, ], nrow=length(keep), ncol=dimension )
        weights <- weights[ keep ]
    }
    
    # If symmetric rule: so far, it's for the positive orthant. 
	# Copy to other orthants!
    if ( sym ) {
        nr  <- length( weights )
        m   <- n1D[[ 1 ]]  # scalar
        for ( j in 1:dimension ) {
            keep    <- rep( 0, nr )
            numnew  <- 0
            for ( r in 1:nr ) {
                if ( nodes[r,j] != m ) {
                    numnew          <- numnew + 1
                    keep[ numnew ]  <- r
                }
            }
            if ( numnew > 0 ) {
                # add as.matrix statement to keep nodes as a matrix, if it contains only 1 column
                nodes                           <- rbind( nodes, matrix( nodes[ keep[ 1:numnew ], ], nrow=numnew, ncol=dimension ) )
                nodes[ (nr+1):(nr+numnew), j ]  <- 2*m - nodes[ (nr+1):(nr+numnew), j ]
                weights                         <- c( weights, weights[ keep[ 1:numnew ] ] )
                nr                              <- nr + numnew
            }
        }
    
        # 3. final sorting
        nodes.sorted    <- sortrows( nodes, index.return=TRUE )    # sortvec
        nodes           <- nodes.sorted$x
        weights         <- weights[ nodes.sorted$ix ]
    }
    
    # normalize weights to account for rounding errors
    weights <- weights / sum(weights)
    
    return( list( "nodes"     = nodes,
                  "weights"   = weights ) )
}

