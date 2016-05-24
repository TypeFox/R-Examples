fuse.Hst.ls <-
function(Hst.ls1, Hst.ls2) {
    tau <- length(Hst.ls1)
    for(i in 1:tau) {
        Hst.ls1[[i]] <- cbind( Hst.ls1[[i]], Hst.ls2[[i]] )
    }
    return(Hst.ls1)
}
