applystnd.Hst.ls <-
function( Hst0.ls, x ) {
    tau <- length(Hst0.ls)
	sHst0.ls <- list()
	for(i in 1:tau) { 
		sHst0.ls[[i]] <- t( ( t(Hst0.ls[[i]]) - x$h.mean ) / x$h.sd )
	}
	return(sHst0.ls)
}
