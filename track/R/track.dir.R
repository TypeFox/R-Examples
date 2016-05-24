track.dir <- function(pos=1, envir=as.environment(pos), data=FALSE) {
    trackingEnv <- getTrackingEnv(envir)
    ## fileMap <- getFileMapObj(trackingEnv)
    ## unsaved <- getUnsavedObj(trackingEnv)
    dir <- getTrackingDir(trackingEnv)
    if (data)
        dir <- getDataDir(dir)
    return(dir)
}

