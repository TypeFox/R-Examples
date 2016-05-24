track.filename <- function(expr, list=character(0), pos=1, envir=as.environment(pos), suffix=FALSE) {
    if (!missing(expr)) {
        ## evaluate expr if necessary, and convert to list
        qexpr <- substitute(expr)
        if (!is.name(qexpr) && !is.character(qexpr))
            stop("expr must be a quoted or unquoted variable name")
        list <- c(as.character(qexpr), list)
    }
    if (!is.character(list))
        stop("'list' must be a character vector")
    trackingEnv <- getTrackingEnv(envir)
    opt <- track.options(trackingEnv=trackingEnv)
    fileMap <- getFileMapObj(trackingEnv)
    i <- match(list, names(fileMap))
    if (suffix)
        return(paste(fileMap[i], opt$RDataSuffix, sep="."))
    else
        return(fileMap[i])
}

track.datadir <- function(pos=1, envir=as.environment(pos), relative=TRUE) {
    trackingEnv <- getTrackingEnv(envir)
    trackingDir <- getTrackingDir(trackingEnv)
    dataDir <- getDataDir(trackingDir)
    if (relative)
        dataDir <- find.relative.path(getwd(), dataDir)
    return(dataDir)
}
