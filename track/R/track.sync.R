track.sync <- function(pos=1, master=c("auto", "envir", "files"), envir=as.environment(pos), trackingEnv=getTrackingEnv(envir), full=TRUE, dryRun=FALSE, taskEnd=FALSE) {
    ## With master="envir", sync the tracking database to the contents of the R environment
    ## This involves 3 things
    ##   (1) start tracking new untracked variables
    ##   (2) for objects that have disappeared from the environment,
    ##       delete them from the tracking database
    ##   (3) check that all variables are activeBindings (to catch cases where
    ##       {rm(x); assign("x", value)} is performed, which leaves tracked
    ##       variables without an active binding.) This is only done as part "full"
    ##   (4) flush excess objects from memory if too much memory is being used
    ##
    ## With master="files", sync the R environment to the tracking database in the filesystem,
    ## which is the job of track.rescan(), so call that.

    ## Could call untracked() and track.orphaned() to find these, but faster
    ## to do ls() on the environment, and look at the tracking file map,
    ## and work out how to sync them up.

    ## Want to use this function as a top-level callback to automatically
    ## track new objects and removed deleted objects, so want it to be fast.

    ## Do check for untrackable objects (isReservedName())
    if (missing(pos)) {
        n <- environmentName(envir)
        if (n=="R_GlobalEnv")
            pos <- 1
        else
            pos <- match(n, search())
    }
    opt <- track.options(trackingEnv=trackingEnv)
    verbose <- dryRun || opt$debug > 0
    trace <- getOption("track.callbacks.trace", FALSE)
    if (verbose)
        cat("track.sync", if (dryRun) "(dryRun)",
            ": syncing tracked env ", envname(envir), "\n", sep="")
    if (!opt$stealable && !opt$readonly) {
        ## The only db we don't want to check for external modifications is a non-stealable writable one
        if (verbose)
            cat('track.sync: proceeding with writeable db ', envname(envir), '\n', sep='')
    } else {
        if (verbose)
            if (opt$stealable)
                cat('track.sync: seeing if stealable db has changed ', envname(envir), '\n', sep='')
            else
                cat('track.sync: seeing if readonly db has changed ', envname(envir), '\n', sep='')
        ## See if the tracking db has changed
        modTime <- file.info(file.path(getTrackingDir(trackingEnv), paste('.trackingSummary', opt$RDataSuffix, sep='.')))
        oldModTime <- mget(envir=trackingEnv, '.trackingModTime', ifnotfound=list(NULL))[[1]]
        if (!is.null(oldModTime) && (modTime$mtime > oldModTime$mtime || modTime$size != oldModTime$size)) {
            cat('track.sync: DB backing ', envname(envir), '[pos=', pos, '] has been modified; rescanning... ', sep='')
            res <- track.rescan(envir=envir, forgetModified=TRUE, level='low', verbose=TRUE)
        }
    }
    master <- match.arg(master)
    if (master=="auto")
        if (opt$readonly)
            master <- "envir"
        else
            stop("must supply argument master='files' or master='envir' when tracking db is attached with readonly=FALSE")
    if (master=="files")
        return(track.rescan(envir=envir, forgetModified=TRUE, level="low"))
    if (trace==2) {
        cat("[", pos, ":", sep="")
        flush.console()
    }

    ## Get info about the state of things
    autoTrack <- mget(".trackAuto", envir=trackingEnv, ifnotfound=list(list(on=FALSE, last=-1)))[[1]]
    fileMap <- getFileMapObj(trackingEnv)
    all.objs <- ls(envir=envir, all.names=TRUE)
    ## 'untracked' will be the untracked vars that we want to track
    untracked <- setdiff(all.objs, names(fileMap))
    reserved <- isReservedName(untracked)
    ## .trackingEnv will always exist -- don't warn about it
    warn.reserved <- setdiff(untracked[reserved], c(".trackingEnv", ".trackingCreated"))
    if (verbose && length(warn.reserved))
        cat("track.sync: cannot track variables with reserved names: ", paste(warn.reserved, collapse=", "), "\n", sep="")
    untracked <- untracked[!reserved]
    if (trace==2 && length(untracked)) {
        cat("u")
        flush.console()
    }
    activeBindings <- sapply(untracked, bindingIsActive, envir)
    if (verbose && any(activeBindings))
        cat("track.sync: cannot track variables that have active bindings: ", paste(untracked[activeBindings], collapse=", "), "\n", sep="")
    untracked <- untracked[!activeBindings]
    if (length(opt$autoTrackExcludeClass)) {
        if (trace==2) {
            cat("e")
            flush.console()
        }
        hasExcludedClass <- sapply(untracked, function(o) any(is.element(class(get(o, envir=envir, inherits=FALSE)), opt$autoTrackExcludeClass)))
        if (any(hasExcludedClass)) {
            if (verbose)
                cat("track.sync: not tracking variables from excluded classes: ",
                    paste(untracked[hasExcludedClass], collapse=", "), "\n", sep="")
            untracked <- untracked[!hasExcludedClass]
        }
    }
    for (re in opt$autoTrackExcludePattern)
        untracked <- grep(re, untracked, invert=TRUE, value=TRUE)
    deleted <- setdiff(names(fileMap), all.objs)

    ## The only thing to do for a readonly env is
    ## to flush cached objects out of memory.
    ## Well..., perhaps we should also check that
    ## no new variables have been created, and if
    ## they have, warn about them.  But, that takes
    ## time, and this function is called after every
    ## top level task...

    ## Deal with new (untracked) and deleted variables
    if (length(untracked)) {
        if (opt$readonly) {
            warning(length(untracked), " variables created in a readonly tracking env, these will not be written out to the files: ",
                    paste(untracked, collapse=", "))
        } else if (dryRun) {
            cat("track.sync(dryRun): would track ", length(untracked), " untracked variables: ", paste(untracked, collapse=", "), "\n", sep="")
        } else {
            if (verbose > 0)
                cat("track.sync: tracking ", length(untracked), " untracked variables: ", paste(untracked, collapse=", "), "\n", sep="")
            if (trace==2) {
                cat("t")
                flush.console()
            }
            track(list=untracked, envir=envir)
            fileMap <- getFileMapObj(trackingEnv)
        }
    } else {
        if (verbose)
            cat("track.sync: no untracked variables\n")
    }
    if (length(deleted)) {
        if (opt$readonly) {
            warning(length(deleted), " variables deleted from a readonly tracking env, these will not be deleted from the files: ",
                    paste(deleted, collapse=", "))
        } else if (dryRun) {
            cat("track.sync(dryRun): would remove ", length(deleted), " deleted variables: ", paste(deleted, collapse=", "), "\n", sep="")
        } else {
            if (verbose > 0)
                cat("track.sync: removing ", length(deleted), " deleted variables: ", paste(deleted, collapse=", "), "\n", sep="")
            if (trace==2) {
                cat("d")
                flush.console()
            }
            track.remove(list=deleted, envir=envir, force=TRUE)
            fileMap <- getFileMapObj(trackingEnv)
        }
    } else {
        if (verbose)
            cat("track.sync: no deleted variables\n")
    }
    now <- as.numeric(proc.time()[3])
    full.orig <- full
    if (is.na(full)) {
        if (opt$autoTrackFullSyncWait==0) {
            full <- TRUE
        } else if (opt$autoTrackFullSyncWait>0) {
            if (autoTrack$last < 0 || now - autoTrack$last >= opt$autoTrackFullSyncWait)
                full <- TRUE
        }
        if (is.na(full))
            full <- FALSE
    }
    if (!opt$readonly) {
        ## Don't look for changes in a readonly db -- takes too long
        ## (there could be changes, and we could warn about them...)
        retrack <- character(0)
        if (full) {
            if (trace==1 && is.na(full.orig)) {
                cat("track.sync.callback", envname(envir), ": look for vars without active bindings at ", date(), "\n", sep="")
                stime <- proc.time()
            }
            if (trace==2) {
                cat("f")
                flush.console()
            }
            ## Find the vars that look like they are tracked but don't have active bindings
            ## This can be time consuming -- need to call bindingIsActive for each tracked
            ## var.
            tracked <- intersect(names(fileMap), all.objs)
            reserved <- isReservedName(tracked)
            if (verbose && any(reserved))
                cat("track.sync: cannot track variables with reserved names: ", paste(tracked[reserved], collapse=", "), "\n", sep="")
            tracked <- tracked[!reserved]
            if (length(tracked))
                retrack <- tracked[!sapply(tracked, bindingIsActive, envir)]
            if (trace==1 && is.na(full.orig)) {
                cat("track.sync.callback: finished looking for vars without active bindings",
                            " (", paste(round(1000*(proc.time()-stime)[1:3]), c("u", "s", "e"), sep="", collapse=" "), " ms)\n", sep="")
            }
        }
        if (length(retrack))
            for (re in opt$autoTrackExcludePattern)
                retrack <- grep(re, retrack, invert=TRUE, value=TRUE)
        ## Deal with untracked objects in the tracked env.
        ## Need to write these to files, and replace with active bindings.
        if (trace==2 && length(retrack)) {
            cat("r")
            flush.console()
        }
        for (objName in retrack) {
            ## get obj from envir, store in file, create active binding
            objval <- get(objName, envir=envir, inherits=FALSE)
            if (any(is.element(class(objval), opt$autoTrackExcludeClass))) {
                if (verbose)
                    cat("track.sync", if (dryRun) "(dryRun)", ": var is from excluded class, not tracking: ", objName, "\n", sep="")
                next
            }
            if (verbose && !opt$readonly)
                cat("track.sync: retracking var: ", objName, "\n", sep="")
            if (opt$readonly)
                warning("binding for variable ", objName, " was clobbered in a readonly tracking env -- forgetting the new value, restoring the old")
            if (dryRun)
                next
            ## Use setTrackedVar to write the object to disk (or merely cache
            ## it in trackingEnv, depending on settings in opt).
            ## setTrackedVar() will assign it in the trackingEnv -- it currently
            ## exists in 'envir'
            if (!opt$readonly)
                setTrackedVar(objName, objval, trackingEnv, opt)
            if (opt$debug >= 2)
                cat('track.sync: removing', paste(objName, collapse=', '), 'from trackingEnv\n')
            remove(list=objName, envir=envir)
            f <- createBindingClosure(objName, trackingEnv)
            makeActiveBinding(objName, env=envir, fun=f)
        }
        ## Do we need to re-read the fileMap?
        if (length(retrack))
            fileMap <- getFileMapObj(trackingEnv)
    }

    ## Which variables are currently cached?
    ## Can't do this until after have checked for untracked vars,
    ## otherwise won't treat those properly.
    ## This code used to call track.flush(envir=envir, all=TRUE)
    ## but that's slow compared to working out flushVars here and
    ## passing the specific vars to track.flush()
    ##
    ## Record variables to flush from cache in flushVars.
    ## Record variables that need saving to disk in saveVars
    ## Note that vars are actually flushed from cache by a
    ## call to track.flush(), which won't flush vars named
    ## in opt$alwaysCache.  For vars named in opt$alwaysCache,
    ## we do want to write them out to file (if they've changed),
    ## but we don't want to remove them from the tracking env.
    if (taskEnd && opt$cachePolicy=="eotPurge") {
        if (trace==2) {
            cat("p")
            flush.console()
        }
        flushVars <- NULL
        saveVars <- NULL
        unsavedVars <- getUnsavedObj(trackingEnv)
        objSummary <- getObjSummary(trackingEnv, opt=opt)
        if (!is.null(objSummary)) {
            ## which variables are currently cached in memory and are candidate for flushing?
            inmem <- is.element(rownames(objSummary), ls(envir=trackingEnv, all.names=TRUE))
            keep1 <- !is.na(objSummary$cache) & (objSummary$cache %in% c("yes", "fixedyes"))
            flushCand <- inmem & !keep1
            if (!any(flushCand)) {
                flushVars <- character(0)
            } else if (length(opt$cacheKeepFun) && !identical(opt$cacheKeepFun, "none")) {
                if (trace==2) {
                    cat("k")
                    flush.console()
                }
                ## If there is a cacheKeepFun, see what it says...
                keep <- try(do.call(opt$cacheKeepFun, list(objs=objSummary, inmem=flushCand, envname=envname(envir))), silent=TRUE)
                ## Expecting a logical vector matching rows of objSummary.
                ## Be informative about any problems with what it returns because this can be a user-supplied function.
                if (is(keep, "try-error")) {
                    warning("opt$cacheKeepFun on ", envname(envir), " stopped with an error: ", keep)
                    keep <- FALSE
                } else if (!is.atomic(keep)) {
                    warning("opt$cacheKeepFun on ", envname(envir), " returned a ", class(keep), " object; expecting a logical vector")
                    keep <- FALSE
                } else if (length(keep)!=nrow(objSummary)) {
                    warning("opt$cacheKeepFun on ", envname(envir), " returned an object of length ", length(keep),
                            "; expected a logical vector of length ", nrow(objSummary))
                    keep <- FALSE
                } else if (is.numeric(keep)) {
                    if (any(is.na(keep) | (keep != 0 & keep != 1)))
                        warning("opt$cacheKeepFun on ", envname(envir), " returned a numeric vector with values other than 0 or 1; interpreting as > 0 as TRUE")
                    keep <- ifelse(is.na(keep), FALSE, keep > 0)
                } else if (!is.logical(keep)) {
                    warning("opt$cacheKeepFun on ", envname(envir), " returned a ", mode(keep), " vector; expecting logical")
                    keep <- FALSE
                } else if (any(is.na(keep))) {
                    warning("opt$cacheKeepFun on ", envname(envir), " returned NAs for ", sum(is.na(keep)), " object(s); e.g.: ",
                            paste(rownames(objSummary)[head(which(is.na(keep)), 3)], collapse=', '))
                    keep <- ifelse(is.na(keep), FALSE, keep)
                }
                flushVars <- rownames(objSummary)[flushCand & !keep]
                saveVars <- intersect(rownames(objSummary)[flushCand & keep], unsavedVars)
            } else {
                keep <- FALSE
                flushVars <- rownames(objSummary)[flushCand & !keep]
                saveVars <- intersect(rownames(objSummary)[flushCand & keep], unsavedVars)
            }
            if (length(flushVars)) {
                if (length(opt$alwaysCache)) {
                    i <- is.element(flushVars, opt$alwaysCache)
                    if (any(i)) {
                        saveVars <- unique(c(saveVars, intersect(flushVars[i], unsavedVars)))
                        flushVars <- flushVars[!i]
                    }
                }
            }
            ## Add the vars we want to write to disk and also keep in memory because
            ## of the cache flag (from alwaysCacheClass)
            saveVars <- unique(c(saveVars, intersect(rownames(objSummary)[keep1 & inmem], unsavedVars)))
        } else {
            warning(".trackingSummary does not exist in trackingEnv ", envname(trackingEnv))
            flushVars <- ls(envir=trackingEnv, all.names=TRUE)
            flushVars <- flushVars[is.element(flushVars, names(fileMap))]
        }
        if (dryRun) {
            cat("track.sync(dryRun): Would flush", length(flushVars), "vars:",
                paste(flushVars, collapse=", "), "\n")
            cat("track.sync(dryRun): Would save", length(saveVars), "vars:",
                paste(saveVars, collapse=", "), "\n")
        } else {
            if (length(flushVars)) {
                if (verbose)
                    cat("track.sync: flushing ", length(flushVars), " vars with call to track.flush(envir=",
                        envname(envir), ", list=c(", paste("'", flushVars, "'", sep="", collapse=", "), "))\n", sep="")
                if (trace==2) {
                    cat("f")
                    flush.console()
                }
                track.flush(envir=envir, list=flushVars)
            }
            if (length(saveVars)) {
                if (verbose)
                    cat("track.sync: saving ", length(saveVars), " vars with call to track.save(envir=",
                        envname(envir), ", list=c(", paste("'", saveVars, "'", sep="", collapse=", "), "))\n", sep="")
                if (trace==2) {
                    cat("s")
                    flush.console()
                }
                track.save(envir=envir, list=saveVars)
            }
        }
    } else if (!opt$readonly) {
        if (dryRun) {
            cat("track.sync(dryRun): Would save all vars\n")
        } else {
            if (verbose)
                cat("track.sync: calling track.save(envir=", envname(envir), ")\n", sep="")
            if (trace==2) {
                cat("s")
                flush.console()
            }
            track.save(envir=envir, all=TRUE)
        }
    }
    if (!opt$readonly) {
        ##  write out the object summary if necessary
        summaryChanged <- mget(".trackingSummaryChanged", ifnotfound=list(FALSE), envir=trackingEnv)[[1]]
        if (summaryChanged) {
            if (trace==2) {
                cat("S")
                flush.console()
            }
            if (!exists(".trackingSummary", envir=trackingEnv, inherits=FALSE)) {
                warning("no .trackingSummary in trackng env ", envname(trackingEnv))
            } else {
                dir <- getTrackingDir(trackingEnv)
                if (verbose)
                    cat("track.sync: saving .trackingSummary for envir=", envname(envir), " to ", dir, "\n", sep="")
                save.res <- saveObjSummary(trackingEnv, opt=opt, dataDir=getDataDir(dir))
                if (is(save.res, "try-error"))
                    warning("unable to save .trackingSummary to ", dir)
                else
                    assign(".trackingSummaryChanged", FALSE, envir=trackingEnv)
            }
        }
    }
    if (full & !dryRun) {
        ## save the time that we did this full sync
        autoTrack$last <- now
        assign(".trackAuto", autoTrack, envir=trackingEnv)
    }
    if (trace==2) {
        cat("]")
        flush.console()
    }
    return(invisible(list(new=untracked, removed=deleted)))
}
