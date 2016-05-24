subsetsites.Hst.ls <-
function( Hst.ls, xmask ) {
    tau <- length(Hst.ls)
	Hst.ls.out <- list()
    
    for(i in 1:tau) { 
		Hst.ls.out[[i]] <-  Hst.ls[[i]][ xmask, , drop=FALSE]
	}
    return( Hst.ls.out )
}
