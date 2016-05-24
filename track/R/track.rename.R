track.rename <- function(old, new, pos=1, envir=as.environment(pos), clobber=FALSE, verbose=TRUE) {
    if (!is.character(old))
        stop("'old' must be a character vector")
    if (!is.character(new) || length(old)!=length(new))
        stop("'new' must be a character vector the same length as 'old'")
    if (any(duplicated(new)))
        stop("'new' has duplicate names: ", paste(duplicated(new), collpase=", "))
    if (length(old)==0)
        return(invisible(list(old=old, new=new)))
    trackingEnv <- getTrackingEnv(envir)
    opt <- track.options(trackingEnv=trackingEnv)
    if (opt$readonly)
        stop("cannot rename in a readonly tracking env")
    all.objs <- ls(envir=envir, all.names=TRUE)
    fileMap <- getFileMapObj(trackingEnv)
    if (any(isReservedName(new)))
        stop("'new' contains reserved names: ", new[isReservedName(new)])
    if (any(isReservedName(old)))
        stop("'old' contains reserved names: ", old[isReservedName(old)])
    if (!clobber && any(new %in% all.objs))
        stop("clobber=FALSE and 'new' names already exist in 'to': ", paste(intersect(new, all.objs), collapse=", "))
    if (!all(old %in% all.objs))
        stop("some 'old' names don't exist in the tracked env:", paste(setdiff(old, all.objs), collapse=", "))
    # check for db consistency
    if (!all(intersect(old, names(fileMap)) %in% all.objs))
        stop("problem with tracked db: some tracked 'old' listed in the fileMap don't exist in the tracked env (run track.rebuild() to try to repair): ",
             paste(setdiff(intersect(old, names(fileMap)), all.objs), collapse=", "))
    # make sure tracked objects are flushed out to files
    track.flush(envir=envir, list=intersect(union(new, old), names(fileMap)))
    objSmy <- getObjSummary(trackingEnv, opt=opt)
    dir <- getDataDir(getTrackingDir(trackingEnv))
    # what makes this tricky:
    #   (1) overlap between 'new' and 'old' names -- need to make temp
    #       copies (of possibly files and/or objects)
    #   (2) untracked 'old' names: keep these as untracked new names (having
    #       this policy means we don't need to worry about exclude classes)
    #   (3) 'new' names that match the exclude pattern
    #   (4) file names for new objects that might change from an underscore
    #       name to a name matching the object name (e.g., changing from
    #       'X1p.v7.2' to 'x')
    old.isTracked <- is.element(old, names(fileMap))
    new.isTracked <- old.isTracked & !exclude.from.tracking(new, opt=opt)
    # Cache 'old' untracked objects whose name also appears in 'new' and
    # old tracked objects that will become untracked.  Reading in objects
    # that will become untracked allows us to fail early before modifying
    # tracking db if we can't fit all the objects in memory.
    old.isCached <- !old.isTracked & is.element(old, new)
    if (any(old.isCached | (!new.isTracked & old.isTracked))) {
        old.cacheEnv <- new.env()
        on.exit(remove(list=ls(envir=old.cacheEnv, all.names=TRUE), envir=old.cacheEnv), add=TRUE)
        for (objName in old[old.isCached | (!new.isTracked & old.isTracked)])
            assign(objName, envir=old.cacheEnv, get(objName, envir=envir, inherits=FALSE))
    }
    # Make a file cache for 'old' tracked objects whose name also appears in 'new'
    # To be safe, make copies of files rather than renaming (which would delete
    # the original before we've successfully created the copy.)
    # Keep this on the same filesystem (rather than using tempdir()) so that
    # we can use file.rename().
    old.inTempFile <- old.isTracked & is.element(old, new)
    if (any(old.inTempFile)) {
        old.tempDir <- tempfile(pattern="__renamedir", tmpdir=dir)
        dir.create(old.tempDir)
        on.exit(unlink(old.tempDir, recursive=TRUE), add=TRUE)
        if (!dir.exists(old.tempDir))
            stop("failed to create temporary directory '", old.tempDir, "' to use for temp copies of files")
        for (objName in old[old.inTempFile]) {
            oldFile <- file.path(dir, paste(fileMap[objName], opt$RDataSuffix, sep="."))
            tmpFile <- file.path(old.tempDir, paste(fileMap[objName], opt$RDataSuffix, sep="."))
            if (!file.copy(oldFile, tmpFile) || !file.exists(tmpFile))
                stop("failed to make temporary copy of file '", oldFile, "'")
            # confirm that we were able to copy the file ok
            info.from <- file.info(oldFile)
            info.to <- file.info(tmpFile)
            if (info.from$size != info.to$size)
                stop("copied file '", oldFile, "' for obj '", objName, "' to '", tmpFile, "', but size of copy (",
                     info.to$size, " bytes) does not match size of original (", info.from$size, " bytes)")
        }
    }
    # Vars that are simply being clobbered (and not renamed themselves)
    # cause complications with updating the fileMap and objSmy, so just
    # remove them now.
    if (clobber) {
        dispose <- intersect(setdiff(new, old), all.objs)
        if (length(dispose)) {
            track.remove(list=dispose, envir=envir)
            fileMap <- getFileMapObj(trackingEnv)
            objSmy <- getObjSummary(trackingEnv, opt=opt)
        }
    }
    # Where we have names in both old and new, can't update fileMap and
    # objSmy one-at-a-time because we will get duplicate indices with
    # the consequent ambiguity.  For these, we need to update at the
    # end, which means we risk a job half done.  So, the intervening
    # code should be as robust as possible to file-system errors, etc.
    obj.updateAtEnd <- (old %in% new) | (new %in% old)
    renamed.ok <- rep(FALSE, length(old))
    tmpEnv <- new.env()
    for (obj.i in seq(along=old)) {
        oldObjName <- old[obj.i]
        newObjName <- new[obj.i]
        if (verbose)
            cat("Renaming", " '", oldObjName, "' to '", newObjName, "'...\n", sep="")
        fileMap.changed <- FALSE
        objSmy.changed <- FALSE
        fileMap.i <- match(oldObjName, names(fileMap))
        objSmy.i <- match(oldObjName, rownames(objSmy))
        # Three situations:
        #    (1) old obj name is tracked and new obj name is tracked
        #    (2) old obj name is tracked and new obj name is not tracked
        #    (3) old obj name is not tracked and new obj name is not tracked
        if (old.isTracked[obj.i]) {
            if (old.inTempFile[obj.i])
                oldFile <- file.path(old.tempDir, paste(fileMap[fileMap.i], opt$RDataSuffix, sep="."))
            else
                oldFile <- file.path(dir, paste(fileMap[fileMap.i], opt$RDataSuffix, sep="."))
            oldFile.needsRemoving <- TRUE
        }
        if (old.isTracked[obj.i] & new.isTracked[obj.i]) {
            # situation (1)
            old.isSimpleName <- isSimpleName(oldObjName)
            new.isSimpleName <- isSimpleName(newObjName)
            # If both old and new are not simple, can keep using the same underlying file,
            # otherwise need to rename
            if (!old.isSimpleName && !new.isSimpleName)
                oldFile.needsRemoving <- FALSE
            newFileName <- makeObjFileName(newObjName, fileMap)
            newFile <- file.path(dir, paste(newFileName, opt$RDataSuffix, sep="."))
            # Unfortunately, can't just rename the file because the object name has changed,
            # and the object is stored in the file.
            # old code doesn't work: ok <- file.rename(oldFile, newFile)
            res <- try(load(oldFile, envir=tmpEnv), silent=TRUE)
            if (is(res, "try-error")) {
                warning("failed to rename obj '", oldObjName, "' to '", newObjName,
                        "' because load('", oldFile, "') failed: ", res)
                next
            }
            if (length(res)!=1 || res!=oldObjName) {
                warning("failed to rename obj '", oldObjName, "' to '", newObjName,
                        "' because old file '", oldFile, "' did not contain object named '", oldObjName, "'")
                next
            }
            if (newObjName != oldObjName) {
                # perverse sitation if newObjName==oldObjName, but try to handle all gracefully...
                res <- try(assign(newObjName, get(oldObjName, envir=tmpEnv, inherits=FALSE), envir=tmpEnv), silent=TRUE)
                if (is(res, "try-error")) {
                    warning("failed to rename obj '", oldObjName, "' to '", newObjName,
                            "' because could not make copy of object for saving")
                    next
                }
                remove(list=oldObjName, envir=tmpEnv, inherits=FALSE)
            } else {
                oldFile.needsRemoving <- FALSE
            }
            res <- try(save(list=newObjName, file=newFile, envir=tmpEnv,
                            compress=opt$compress, compression_level=opt$compression_level), silent=TRUE)
            if (is(res, "try-error")) {
                warning("failed to rename obj '", oldObjName, "' to '", newObjName,
                        "' because could not save renamed obj in file '", oldFile, "': ", res)
                next
            }
            remove(list=newObjName, envir=tmpEnv, inherits=FALSE)
            fileMap[fileMap.i] <- newFileName
            if (!obj.updateAtEnd[obj.i]) {
                names(fileMap)[fileMap.i] <- newObjName
                fileMap.changed <- TRUE
                if (!is.na(objSmy.i)) {
                    rownames(objSmy)[objSmy.i] <- newObjName
                    objSmy.changed <- TRUE
                }
            }
            # remove any cached object in the 'to' envir
            if (exists(newObjName, envir=trackingEnv, inherits=FALSE))
                remove(list=newObjName, envir=trackingEnv)
            # and remove the active binding if it exists
            if (exists(newObjName, envir=envir, inherits=FALSE))
                remove(list=newObjName, envir=envir)
            f <- createBindingClosure(newObjName, trackingEnv)
            makeActiveBinding(newObjName, env=envir, fun=f)
        } else if (!new.isTracked[obj.i]) {
            # situations (2) & (3) (new obj is not tracked)
            if (old.isCached[obj.i])
                objVal <- get(oldObjName, envir=old.cacheEnv, inherits=FALSE)
            else
                objVal <- get(oldObjName, envir=envir, inherits=FALSE)
            assign(newObjName, envir=envir, value=objVal)
            if (!is.na(objSmy.i)) {
                objSmy <- objSmy[ -objSmy.i, , drop=FALSE]
                objSmy.changed <- TRUE
            }
        } else {
            stop("shouldn't happen: old objName is not tracked and new objName is tracked: old='",
                 oldObjName, "', new='", newObjName, "'")
        }
        if (objSmy.changed) {
            assign.res <- try(assign(".trackingSummary", objSmy, envir=trackingEnv), silent=TRUE)
            if (is(assign.res, "try-error")) {
                stop("unable to assign .trackingSummary back to tracking env on ",
                        envname(trackingEnv), ": ", assign.res)
            } else {
                assign(".trackingSummaryChanged", TRUE, envir=trackingEnv)
                save.res <- saveObjSummary(trackingEnv, opt=opt, dataDir=getDataDir(dir))
                if (is(save.res, "try-error"))
                    stop("unable to save .trackingSummary to ", dir)
                else
                    assign(".trackingSummaryChanged", FALSE, envir=trackingEnv)
            }
        }
        if (fileMap.changed) {
            writeFileMapFile(fileMap, trackingEnv=trackingEnv, dataDir=getDataDir(dir))
        }
        if (!oldObjName %in% new)
            remove(list=oldObjName, envir=envir, inherits=FALSE)
        if (old.isTracked[obj.i] && oldFile.needsRemoving) {
            if (!file.remove(oldFile))
                warning("unable to remove file '", oldFile, "'")
        }
        renamed.ok[obj.i] <- TRUE
    }
    if (any(obj.updateAtEnd)) {
        # these are the ones we couldn't do one at a time because they were in both
        # old and new
        names(fileMap)[match(old[obj.updateAtEnd], names(fileMap))] <- new[obj.updateAtEnd]
        rownames(objSmy)[match(old[obj.updateAtEnd], rownames(objSmy))] <- new[obj.updateAtEnd]
        assign.res <- try(assign(".trackingSummary", objSmy, envir=trackingEnv), silent=TRUE)
        if (is(assign.res, "try-error")) {
            stop("unable to assign .trackingSummary back to tracking env on ",
                    envname(trackingEnv), ": ", assign.res)
        } else {
            assign(".trackingSummaryChanged", TRUE, envir=trackingEnv)
            save.res <- saveObjSummary(trackingEnv, opt=opt, dataDir=getDataDir(dir))
            if (is(save.res, "try-error"))
                stop("unable to save .trackingSummary to ", dir)
            else
                assign(".trackingSummaryChanged", FALSE, envir=trackingEnv)
        }
        writeFileMapFile(fileMap, trackingEnv=trackingEnv, dataDir=getDataDir(dir))
    }
    return(invisible(list(old=old, new=new)))
}
