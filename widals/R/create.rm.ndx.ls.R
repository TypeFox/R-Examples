create.rm.ndx.ls <-
function(n, xincmnt=10) {
    rm.ndx.ls <- list()
    for(i in 1:xincmnt) { xrm.ndxs <- seq(i, n+xincmnt, by=xincmnt) ; xrm.ndxs <- xrm.ndxs[ xrm.ndxs <= n ] ; rm.ndx.ls[[i]] <- xrm.ndxs }
    return( rm.ndx.ls )
}
