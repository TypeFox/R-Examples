track.load <- function(files, pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, cache=FALSE, clobber=FALSE, time.of.file=TRUE, warn=TRUE) {
    ## load some or all variables from some RData files into a tracked environment
    trackingEnv <- getTrackingEnv(envir)
    opt <- track.options(trackingEnv=trackingEnv)
    if (!is.logical(cache) || is.na(cache))
        stop("'cache' must be TRUE or FALSE")
    opt$cache <- cache
    tmpenv <- new.env(hash=TRUE, parent=emptyenv())
    all.loaded <- character(0)
    all.skipped <- character(0)
    all.wanted <- list
    for (file in files) {
        if (!file.exists(file)) {
            warning("skipping ", file, ": file does not exist")
            next
        }
        load.res <- try(load(file=file, envir=tmpenv), silent=TRUE)
        if (is(load.res, "try-error")) {
            warning("skipping ", file, ": could not read saved objects: ", as.character(load.res))
            next
        }
        info <- if (time.of.file) file.info(file) else NULL
        list <- all.wanted
        if (is.null(list))
            list <- load.res
        else
            list <- intersect(load.res, list)
        i <- isReservedName(list)
        if (any(i))
            list <- list[!i]
        if (!is.null(glob))
            pattern <- glob2rx(glob)
        if (!is.null(pattern))
            list <- grep(pattern, list, value=TRUE)
        all.objs <- ls(envir=envir, all.names=TRUE)
        i <- is.element(list, all.objs) | duplicated(list)
        if (any(i)) {
            if (!clobber) {
                warning("skipping ", sum(i), " variable(s) because these exist in ", envname(envir),
                        " and clobber=FALSE: ",
                        paste("'", list[i][seq(len=min(3, sum(i)))], "'", sep="", collapse=", "), if (sum(i)>3) "...", "\n")
                list <- list[!i]
            } else {
                warning("clobbering ", sum(i), " existing variable(s) in ", envname(envir),
                        ": ",
                        paste("'", list[i][seq(len=min(3, sum(i)))], "'", sep="", collapse=", "), if (sum(i)>3) "...", "\n")
            }
        }
        for (objName in list) {
            objval <- get(objName, envir=tmpenv)
            setTrackedVar(objName, objval, trackingEnv, opt, info)
            ## always remove the object and reassign the binding
            if (is.element(objName, all.objs))
                remove(list=objName, envir=envir)
            f <- createBindingClosure(objName, trackingEnv)
            makeActiveBinding(objName, env=envir, fun=f)
        }
        all.loaded <- c(all.loaded, list)
        all.skipped <- c(all.skipped, setdiff(load.res, list))
    }
    if (any(!is.element(all.wanted, all.loaded)))
        warning("the following requested objects were not found: ", paste(setdiff(all.wanted, all.loaded), collapse=", "))

    return(list(loaded=all.loaded, skipped=all.skipped))
}

