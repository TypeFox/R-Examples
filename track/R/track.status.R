tracked <- function(pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, all.names=TRUE)
    track.status(envir=envir, list=list, pattern=pattern, glob=glob, file.status=FALSE, what="tracked", all.names=all.names)

untracked <- function(pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, all.names=TRUE)
    track.status(envir=envir, list=list, pattern=pattern, glob=glob, file.status=FALSE, what="untracked", tracked=FALSE, all.names=all.names)

track.orphaned <- function(pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, all.names=TRUE)
    track.status(envir=envir, list=list, pattern=pattern, glob=glob, file.status=FALSE, what="orphaned", tracked=TRUE, all.names=all.names)

track.masked <- function(pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, all.names=TRUE)
    track.status(envir=envir, list=list, pattern=pattern, glob=glob, file.status=FALSE, what="masked", tracked=TRUE, all.names=all.names)

untrackable <- function(pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, all.names=TRUE)
    track.status(envir=envir, list=list, pattern=pattern, glob=glob, file.status=FALSE, what="untrackable", tracked=FALSE, reserved=TRUE, all.names=all.names)

track.unsaved <- function(pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, all.names=TRUE)
    track.status(envir=envir, list=list, pattern=pattern, glob=glob, file.status=FALSE, what="unsaved", tracked=TRUE, all.names=all.names)

track.status <- function(pos=1, envir=as.environment(pos), expr, qexpr=NULL, list=NULL,
                         pattern=NULL, glob=NULL, file.status=TRUE,
                         tracked=NA, reserved=FALSE, all.names=FALSE,
                         what=c("all", "tracked", "trackable", "untracked", "orphaned", "masked", "unsaved", "untrackable")) {
    ## return a dataframe with the status on each var matching pattern or glob, with
    ## the following columns:
    ##   rowname: object name
    ##   status: "tracked", "untracked", "untrackable", "orphaned", "masked", "nonexistent"
    ##   fileBase: base name of file, if status is "tracked"
    ##   inMem: TRUE/FALSE depending on whether obj exists in the tracking env, NA if not tracked
    ##   saved: TRUE/FALSE depending on whether the version in memory has been saved to disk
    ##
    ## status:
    ##   tracked: obj exists in the fileMap and is an active binding on envir
    ##   orphaned: obj exists in the fileMap but does not have any binding on envir
    ##   masked: obj exists in the fileMap but has an ordinary binding on envir
    ##   untracked: obj does not exists in the fileMap and is not an active binding on envir
    ##   untrackable: obj does not exists in the fileMap and is an active binding on envir (and hence
    ##     cannot be tracked)
    ##   nonexistent: obj does not exist in the fileMap and does not exist in envir
    what <- match.arg(what)
    trackingEnv <- getTrackingEnv(envir)
    opt <- track.options(trackingEnv=trackingEnv)
    fileMap <- getFileMapObj(trackingEnv)
    unsaved <- getUnsavedObj(trackingEnv)
    dir <- getTrackingDir(trackingEnv)
    dataDir <- getDataDir(dir)

    if (!missing(expr))
        qexpr <- substitute(expr)
    if (!is.null(qexpr)) {
        if (is.name(qexpr) || is.character(qexpr)) {
            objName <- as.character(qexpr)
        } else {
            stop("expr argument must be a quoted or unquoted variable name")
        }
        list <- c(objName, list)
    }
    all.objs <- ls(envir=envir, all.names= all.names || reserved)
    fileMap.names <- names(fileMap)
    if (!(all.names || reserved))
        fileMap.names <- fileMap.names[substring(fileMap.names, 1, 1)!="."]
    if (is.null(list)) {
        list <- union(all.objs, fileMap.names)
        if (!reserved) {
            i <- isReservedName(list)
            if (any(i))
                list <- list[!i]
        }
        if (!is.null(glob))
            pattern <- glob2rx(glob)
        if (!is.null(pattern))
            list <- grep(pattern, list, value=TRUE)
    } else {
        if (!is.null(glob) || !is.null(pattern))
            stop("must specify EITHER expr or list OR glob or pattern")
    }
    if (!is.na(tracked))
        if (tracked)
            list <- intersect(list, fileMap.names)
        else
            list <- setdiff(list, fileMap.names)
    status <- lapply(list, function(objName) {
        ## an object is already tracked if the following 2 conditions are met:
        ##   - it exists as an activing binding in envir
        ##   - there is an entry in the fileMap in the trackingEnv
        inMem <- fileExists <- as.logical(NA)
        fileBase <- as.character(NA)
        if (is.element(objName, fileMap.names)) {
            saved <- !is.element(objName, unsaved)
            if (!is.element(objName, all.objs))
                status <- "orphaned"
            else if (!bindingIsActive(objName, envir))
                status <- "masked"
            else
                status <- "tracked"
            fileBase <- fileMap[objName]
            if (file.status)
                fileExists <- file.exists(file.path(dataDir, paste(fileBase, opt$RDataSuffix, sep=".")))
            inMem <- exists(objName, trackingEnv, inherits=FALSE)
        } else {
            if (!is.element(objName, all.objs))
                status <- "nonexistent"
            else if (!bindingIsActive(objName, envir) && !isReservedName(objName))
                status <- "untracked"
            else
                status <- "untrackable"
            saved <- as.logical(NA)
        }
        return(list(status, inMem, fileBase, fileExists, saved))
    })
    if (what=="all") {
        if (length(status)) {
            status <- as.data.frame(lapply(1:5, function(i) sapply(status, "[[", i)))
            rownames(status) <- list
        } else {
            status <- data.frame(character(0), logical(0), character(0), logical(0), logical(0))
        }
        colnames(status) <- c("Status", "InMem", "FileBase", "FileExists", "Saved")
        attr(status, "trackingDir") <- dir
        status <- status[order(rownames(status)), , drop=FALSE]
        return(status)
    } else if (what=="unsaved") {
        if (length(list)) {
            saved <- sapply(status, "[[", 5)
            return(sort(list[!saved & !is.na(saved)]))
        } else {
            return(character(0))
        }
    } else {
        if (length(list)) {
            status <- sapply(status, "[[", 1)
            return(sort(list[status==what]))
        } else {
            return(character(0))
        }
    }
}
