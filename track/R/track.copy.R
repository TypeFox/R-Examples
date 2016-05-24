track.copy <- function(from, to=1, list=NULL, pattern=NULL,
                       glob=NULL, delete=FALSE, clobber=FALSE,
                       skipExisting=FALSE,
                       verbose=TRUE, do.untrackable=FALSE) {
    if (!is.numeric(from) && !is.character(from))
        stop("only implemented for numeric or char values for 'from'")
    if (!is.numeric(to) && !is.character(to))
        stop("only implemented for numeric or char values for 'to'")
    if (length(from) > 1) {
        res <- vector('list', length(from))
        names(res) <- as.character(from)
        for (i in seq(along=from)) {
            x <- try(track.copy(from=from[[i]], to=to, list=list, pattern=pattern, glob=glob, delete=delete,
                                clobber=clobber, skipExisting=skipExisting, verbose=verbose, do.untrackable=do.untrackable))
            if (!is(x, 'try-error'))
                res[[i]] <- x
        }
        return(invisible(NULL))
    }
    if (do.untrackable)
        warning("do.untrackable not yet implemented")
    env.to <- as.environment(to)
    env.from <- as.environment(from)
    if (identical(env.to, env.from))
        stop("'from' and 'to' are the same")
    trackingEnv.to <- getTrackingEnv(env.to)
    trackingEnv.from <- getTrackingEnv(env.from)
    if (verbose)
        cat(if (delete) "Moving" else "Copying", " objects from ",
            track.datadir(env.from), " to ",
            track.datadir(env.to), "\n", sep="")
    opt.to <- track.options(trackingEnv=trackingEnv.to)
    opt.from <- track.options(trackingEnv=trackingEnv.from)
    if (opt.to$readonly)
        stop("cannot copy into a readonly tracking env '", to, "'")
    if (opt.from$readonly && delete)
        stop("cannot move (i.e., copy & delete) from readonly tracking env '", from, "'")
    all.objs.from <- ls(envir=env.from, all.names=TRUE)
    all.objs.from <- all.objs.from[!isReservedName(all.objs.from)]
    all.objs.from <- setdiff(all.objs.from, c(".Last", ".Last.sys"))
    fileMap.from <- getFileMapObj(trackingEnv.from)
    # don't limit to just tracked objs
    # all.objs.from <- intersect(all.objs.from, names(fileMap.from))
    if (is.null(pattern) && is.null(glob) && is.null(list))
        list <- all.objs.from
    if (!is.null(pattern))
        list <- unique(c(list, grep(pattern, all.objs.from, value=TRUE)))
    if (!is.null(glob))
        list <- unique(c(list, grep(glob2rx(glob), all.objs.from, value=TRUE)))
    if (!is.null(list)) {
        if (any(isReservedName(list))) {
            warning("omitting objects with reserved names: ", paste(list[isReservedName(list)], collapse=", "))
            list <- list[!isReservedName(list)]
        }
        if (any(!(list %in% all.objs.from))) {
            warning("omitting objects not present in 'from': ", paste(setdiff(list, all.objs.from), collapse=", "))
            list <- setdiff(list, all.objs.from)
        }
    }
    if (length(list)==0) {
        if (verbose)
            cat("Nothing to copy\n")
        return(invisible(list))
    }
    fileMap.to <- getFileMapObj(trackingEnv.to)
    all.objs.to <- ls(envir=env.to, all.names=TRUE)
    if (skipExisting)
        list <- setdiff(list, all.objs.to)
    if (!clobber && any(list %in% all.objs.to))
        stop("clobber=FALSE and some objects to be copied already exist in 'to': ", paste(intersect(list, all.objs.to), collapse=", "))
    # make sure objects in the source are flushed out to files
    track.flush(envir=env.from, list=intersect(list, names(fileMap.from)))
    objSmy.to <- getObjSummary(trackingEnv.to, opt=opt.to)
    objSmy.from <- getObjSummary(trackingEnv.from, opt=opt.from)
    dir.to <- getDataDir(getTrackingDir(trackingEnv.to))
    dir.from <- getDataDir(getTrackingDir(trackingEnv.from))
    auto.to <- track.auto(envir=env.to)
    for (objName in list) {
        if (verbose)
            cat(if (delete) "Moving" else "Copying", " '", objName, "'...\n", sep="")
        fileMap.to.changed <- FALSE
        fileMap.from.i <- match(objName, names(fileMap.from))
        trackedInSource <- !is.na(fileMap.from.i)
        if (trackedInSource) {
            file.from <- fileMap.from[fileMap.from.i]
            file.from.abs <- file.path(dir.from, paste(file.from, opt.from$RDataSuffix, sep="."))
            objSmy.from.i <- match(objName, rownames(objSmy.from))
            objClasses <- strsplit(objSmy.from[objSmy.from.i, "class"], ",")[[1]]
            smyRow <- objSmy.from[objSmy.from.i, ]
        } else {
            objSmy.from.i <- NA
            objVal <- get(objName, envir=env.from, inherits=FALSE)
            objClasses <- class(objVal)
            smyRow <- summaryRow(objName, opt.from, obj=objVal, file=NULL, change=TRUE, times=NULL)
        }
        trackInDest <- (trackedInSource || auto.to) && !exclude.from.tracking(objName, objClasses, opt.to)
        if (trackInDest) {
            if (is.na(file.to <- fileMap.to[match(objName, names(fileMap.to))])) {
                file.to <- makeObjFileName(objName, fileMap.to)
                fileMap.to[objName] <- file.to
                fileMap.to.Changed <- TRUE
            }
            file.to.abs <- file.path(dir.to, paste(file.to, opt.to$RDataSuffix, sep="."))
            if (file.exists(file.to.abs) && !clobber)
                stop("file '", file.to.abs, "' for obj '", objName, "' already exists and clobber=FALSE")
        } else {
            objVal <- get(objName, envir=env.from, inherits=FALSE)
        }
        deleteThis <- FALSE
        if (trackedInSource && trackInDest) {
            # move/copy the file
            if (delete) {
                deleteThis <- TRUE
                # rename is nicest, but won't work across file systems, so just try it
                ok <- file.rename(file.from.abs, file.to.abs)
                if (!ok) {
                    ok <- file.copy(file.from.abs, file.to.abs, overwrite=TRUE)
                    if (!ok)
                        stop("could not copy file '", file.to.abs, "' for obj '", objName, "'")
                    # confirm that we were able to copy the file ok
                    if (!file.exists(file.to.abs))
                        stop("failed to copy file '", file.to.abs, "' for obj '", objName, "'")
                    info.from <- file.info(file.from.abs)
                    info.to <- file.info(file.to.abs)
                    if (info.from$size != info.to$size)
                        stop("copied file '", file.to.abs, "' for obj '", objName, "', but size of copy (",
                             info.to$size, " bytes) does not match size of original (", info.from$size, " bytes)")
                    ok <- file.remove(file.from.abs)
                    if (!ok) {
                        warning("could not remove file '", file.from.abs, "' for obj '", objName, "'")
                        deleteThis <- FALSE
                    }
                }
            } else {
                ok <- file.copy(file.from.abs, file.to.abs, overwrite=TRUE)
                if (!ok)
                    stop("could not copy file '", file.to.abs, "' for obj '", objName, "'")
            }
        } else {
            deleteThis <- delete
        }
        if (trackInDest) {
            # remove any cached object in the 'to' envir
            if (exists(objName, envir=trackingEnv.to, inherits=FALSE))
                remove(list=objName, envir=trackingEnv.to)
            # and remove the active binding if it exists
            if (exists(objName, envir=env.to, inherits=FALSE))
                remove(list=objName, envir=env.to)
            f <- createBindingClosure(objName, trackingEnv.to)
            makeActiveBinding(objName, env=env.to, fun=f)
            if (!trackedInSource) {
                # modify options we give here to immediately write out to disk
                setTrackedVar(objName, value=objVal, trackingEnv=trackingEnv.to, file=file.to,
                              opt=replace(opt.to, c("writeToDisk", "cachePolicy", "maintainSummary"),
                                          list(TRUE, 'none', FALSE)))
            }
        }
        if (!trackInDest)
            assign(objName, objVal, envir=env.to)
        if (trackInDest) {
            # update & write objSmy.to & write fileMap.to
            if (objName %in% rownames(objSmy.to)) {
                objSmy.to[objName, ] <- smyRow
            } else {
                objSmy.to <- rbind(objSmy.to, smyRow)
            }
            assign.res <- try(assign(".trackingSummary", objSmy.to, envir=trackingEnv.to), silent=TRUE)
            if (is(assign.res, "try-error")) {
                stop("unable to assign .trackingSummary back to tracking env on ",
                        envname(trackingEnv.to), ": ", assign.res)
            } else {
                assign(".trackingSummaryChanged", TRUE, envir=trackingEnv.to)
                save.res <- saveObjSummary(trackingEnv=trackingEnv.to, opt=opt.to, dataDir=getDataDir(dir.to))
                if (is(save.res, "try-error"))
                    stop("unable to save .trackingSummary to ", dir.to)
                else
                    assign(".trackingSummaryChanged", FALSE, envir=trackingEnv.to)
            }
            writeFileMapFile(fileMap.to, trackingEnv=trackingEnv.to, dataDir=getDataDir(dir.to))
        }
        if (trackedInSource && deleteThis) {
            # update fileMap.from & objSmy.from
            # can only get here if opt.from$readonly==TRUE
            if (is.na(fileMap.from.i)) {
                warning("tracked object '", objName, "' doesn't have fileMap entry??? - not deleting")
            } else {
                objSmy.from <- objSmy.from[-objSmy.from.i, , drop=FALSE]
                assign.res <- try(assign(".trackingSummary", objSmy.from, envir=trackingEnv.from), silent=TRUE)
                if (is(assign.res, "try-error")) {
                    stop("unable to assign .trackingSummary back to tracking env on ",
                            envname(trackingEnv.from), ": ", assign.res)
                } else {
                    assign(".trackingSummaryChanged", TRUE, envir=trackingEnv.from)
                    save.res <- saveObjSummary(trackingEnv=trackingEnv.from, opt=opt.from, dataDir=getDataDir(dir.from))
                    if (is(save.res, "try-error"))
                        stop("unable to save .trackingSummary to ", dir.from)
                    else
                        assign(".trackingSummaryChanged", FALSE, envir=trackingEnv.from)
                }
                fileMap.from <- fileMap.from[-fileMap.from.i]
                writeFileMapFile(fileMap.from, trackingEnv=trackingEnv.from, dataDir=getDataDir(dir.from))
                if (exists(objName, envir=env.from, inherits=FALSE))
                    remove(list=objName, envir=env.from)
                if (exists(objName, envir=trackingEnv.from, inherits=FALSE))
                    remove(list=objName, envir=trackingEnv.from)
            }
        }
        if (!trackedInSource && deleteThis) {
            remove(list=objName, envir=env.from, inherits=FALSE)
        }
    }
    return(invisible(list))
}

track.move <- function(from, to=1, list=NULL, pattern=NULL, glob=NULL, delete=TRUE, clobber=FALSE, skipExisting=FALSE, verbose=TRUE, do.untrackable=FALSE)
    track.copy(from=from, to=to, list=list, pattern=pattern, glob=glob, delete=delete, clobber=clobber, skipExisting=skipExisting, verbose=verbose, do.untrackable=do.untrackable)
