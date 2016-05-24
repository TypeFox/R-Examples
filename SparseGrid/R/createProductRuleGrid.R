# Copyright (C) 2011 Jelmer Ypma. All Rights Reserved.
# This code is published under the GPL.
#
# File:   createProductRuleGrid.R
# Author: Jelmer Ypma
# Date:   17 May 2011
#
# Input: 
# Output: 
#
# 14/05/2012: Fixed a bug related to checking of arguments when 'type' is a user-supplied function.
#
# TODO: check what happens when we have 1 dimension, should not return anything?
# TODO: works only on built-in functions
# TODO: check for sum(weights)==1 and length(weights)==nrow(nodes)

# create integration grid according to product rule
createProductRuleGrid <- function( type, dimension, k, sym = FALSE ) {
    
    # assertions
    if (! as.integer(dimension) == dimension) {
		stop( paste("createProductRuleGrid expects dimension to be integer, user supplied dimension = ", dimension, sep='') )
	}
	if (! as.integer(k) == k) {
		stop( paste("createProductRuleGrid expects k to be integer, user supplied k = ", k, sep='') )
	}
	if (! is.logical(sym) | is.na(sym)) {
		stop( paste("createProductRuleGrid expects sym to be logical, user supplied sym = ", sym, sep='') )
	}
    
	# check if a function or a string was supplied
    # builtinfct is a boolean defining whether 'type' is a built-in function (TRUE) or a user-supplied function (FALSE)
    builtinfct <- FALSE
	if ( is.function(type) ) {
		if (length( formals( type ) ) != 1 ) {
			stop( "User supplied function (argument type) to createProductRuleGrid needs to have one argument." )
		}
		tmp.grid <- type( dimension )
		if (! is.list(tmp.grid) ) {
			stop( "User supplied function (argument type) to createProductRuleGrid needs to return a list." )
		}
		if (! is.list(tmp.grid) ) {
			stop( "User supplied function (argument type) to createProductRuleGrid needs to return a list." )
		}
		if (! "nodes" %in% names(tmp.grid) ) {
			stop( "Element 'nodes' not found in list returned by user supplied function (argument type) to createProductRuleGrid." )	
		}
		if (! "weights" %in% names(tmp.grid) ) {
			stop( "Element 'weights' not found in list returned by user supplied function (argument type) to createProductRuleGrid." )	
		}
	}
    else if ( is.character(type) ) {
        # input string is from the list of built-in functions?
        builtinfct <- type %in% c('GQU', 'GQN', 'KPU', 'KPN' )
        if ( builtinfct ) { 
            sym = TRUE
        } else {
            stop( paste("createProductRuleGrid expects type to be a string with values 'GQU', 'GQN', 'KPU', 'KPN', or a function, user supplied type = ", type, sep='') )
        }
    }
    else {
        stop( paste("createProductRuleGrid expects type to be a string with values 'GQU', 'GQN', 'KPU', 'KPN', or a function, user supplied type = ", type, sep='') )
    }
    
    grid1d <- createSparseGrid( type, 1, k )

    num.nodes <- length(grid1d$nodes)

    tmp.grid <- grid1d

	if ( num.nodes > 1 ) {
		if ( dimension > 1 ) {
			for (i in 2:dimension) {
				tmp.grid$nodes <- cbind( tmp.grid$nodes[ rep(1:nrow(tmp.grid$nodes), each=num.nodes ), ], rep( grid1d$nodes, times=num.nodes ) )
				tmp.grid$weights <- rep( tmp.grid$weights, each=num.nodes ) * rep( grid1d$weights, times=num.nodes )
			}
		}
		tmp.grid$weights <- tmp.grid$weights / sum( tmp.grid$weights )
	}
	else {
		tmp.grid$nodes <- matrix( rep( grid1d$nodes, dimension ), nrow=1, ncol=dimension )
		tmp.grid$weights <- 1
	}
    
    
    return( tmp.grid )
}
