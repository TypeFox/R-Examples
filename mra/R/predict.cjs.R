predict.cjs <- function( object, ... ){

if( !is.null(object$fitted) ){
	return( object$fitted )
}

f.active <- function( hist ){
	# given a single capture history, return FALSE in the non-active cells,
	# TRUE in the active cells. 
	# In hist, 0 = not captured, 1 = captured and released alive, 2 = 
	# captured and died.  A '2' preceeding a '1' is an error and the latter 1 is ignored 
	# (i.e., the history is stopped at the 2). 
	ns <- length(hist)
	active <- rep( FALSE, ns )
    if( any( hist >= 1 ) ){   # had to add this because all zero rows are allowed
    	first <- min( (1:ns)[hist >= 1] )
    	if( first < ns ){
    		last <- min( ns, (1:ns)[hist == 2] )   # if there are 2 2's, use the first one. 
    		active[ (first+1):last ] <- TRUE
    	}
    }
	active
}
f.cumprod <- function( x ){
	# return cummulative product of elements in x, ignoring NA's
	na.ind <- is.na(x)
	ans <- x
	ans[ !na.ind ] <- cumprod( x[!na.ind] )
	ans
}
act.cells <- t(apply( object$histories, 1, f.active ))
s.hat <- cbind( NA, object$s.hat[,-ncol(object$s.hat)] ) # Put the column of NA's on the front, not back of s.hat
s.hat[ !act.cells ] <- NA 
p.surv.to <- t(apply( s.hat, 1, f.cumprod ))
fitted <- p.surv.to * object$p.hat

fitted
}


