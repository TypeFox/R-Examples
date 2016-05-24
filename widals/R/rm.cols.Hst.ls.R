rm.cols.Hst.ls <-
function(Hst.ls, rm.col.ndx) {
    tau <- length(Hst.ls)
    for(i in 1:tau) {
        Hst.ls[[i]] <- Hst.ls[[i]][  , -rm.col.ndx, drop=FALSE ]
    }
    return(Hst.ls)
}
