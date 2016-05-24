combiner <- function(...) {

    objs <- list(...)
    n <- length(objs)

    if(n == 1) {

	objs <- objs[[ 1 ]]
	n <- length(objs)

    }

    cl <- unlist(lapply(objs, class))
    if(!all(is.element(cl, c("features", "matched")))) stop("combiner: one or more invalid arguments.")

    out <- objs[[ 1 ]]

    if(n > 1) {

	times <- c(unlist(lapply(objs, function(x) attributes(x)$time.point))) 
        attr(out, "time.point") <- times

	xdim <- dim(out$X)

	X.feats <- out$X.feats
	Y.feats <- out$Y.feats

	Xlab <- out$X.labeled
	Ylab <- out$Y.labeled

	X <- out$X
	Xhat <- out$Xhat

	for(i in 2:n) {

	    X.feats <- c(X.feats, objs[[ i ]][[ "X.feats" ]])
	    Y.feats <- c(Y.feats, objs[[ i ]][[ "Y.feats" ]])

	    X <- array(c(c(X), c(objs[[ i ]][[ "X" ]])), dim = c(xdim, i + 1))
	    Xhat <- array(c(c(Xhat), c(objs[[ i ]][[ "Xhat" ]])), dim = c(xdim, i + 1))

	    Xlab <- array(c(c(Xlab), c(objs[[ i ]][[ "X.labeled" ]])), dim = c(xdim, i + 1))
            Ylab <- array(c(c(Ylab), c(objs[[ i ]][[ "Y.labeled" ]])), dim = c(xdim, i + 1))

	} # end of for 'i' loop.

	out$X <- X
	out$Xhat <- Xhat

	out$X.feats <- X.feats
	out$Y.feats <- Y.feats

	out$X.labeled <- Xlab
	out$Y.labeled <- Ylab

    } else {

	warning("combiner: nothing to combine.  Returning input as is.")
	return(out)

    }

    class(out) <- "combined"
    return(out)

} # end of 'combiner' function.
