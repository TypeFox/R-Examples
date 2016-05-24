track.stop <- function(pos=1, envir=as.environment(pos), all=FALSE, stop.on.error=FALSE, keepVars=FALSE, sessionEnd=FALSE, verbose=TRUE, detach=TRUE, callFrom=NULL) {
    force.detach <- FALSE
    if (is.character(detach) && detach=="force") {
        force.detach <- TRUE
        detach <- TRUE
    }
    debug <- !identical(getOption("track.debug", sessionEnd), FALSE)
    ## track.stop() with no arguments behaves analogously
    ## to track.start() with no args, and works on pos=1 (globalenv)
    ## if (missing(pos) && missing(envir) && missing(all)) {
    ##    stop("must specify one of pos, envir, or all")
    ## Want the "sessionEnd" arg for when this is run at the end of an R session --
    ## no recovery possible in that case.
    if (keepVars && sessionEnd) {
        warning("cannot keepVars when sessionEnd==TRUE")
        keepVars <- FALSE
    }
    if (stop.on.error && sessionEnd) {
        warning("cannot use stop.on.error when sessionEnd==TRUE")
        stop.on.error <- FALSE
    }
    ## Detach a tracking env -- should call track.flush first.
    if (all) {
        env.names <- tracked.envs()
        if (verbose && length(env.names) && !is.null(callFrom))
            cat("Stopping all tracking via call from ", callFrom, "\n", sep="")
        for (e.name in env.names)
            if (stop.on.error)
                track.stop(envir=as.environment(e.name), keepVars=keepVars, sessionEnd=sessionEnd)
            else
                try(track.stop(envir=as.environment(e.name), keepVars=keepVars, sessionEnd=sessionEnd, stop.on.error=FALSE))
    } else {
        if (!env.is.tracked(envir)) {
            if (sessionEnd)
                warning("strange: environment ", environmentName(envir), " is not tracked")
            else
                warning("environment ", environmentName(envir), " is not tracked")
            return(invisible(NULL))
        }
        if (verbose)
            cat("Stopping tracking on ", envname(envir), "\n", sep="")
        trackingEnv <- getTrackingEnv(envir)
        opt <- track.options(trackingEnv=trackingEnv)
        if (opt$readonly && keepVars && environmentIsLocked(envir))
            stop("cannot keep vars in a locked readonly tracking environment")
        ## Remove the active bindings, and replace with ordinary variable if keepVars=TRUE
        tracked.vars <- tracked(envir=envir)
        if (!opt$readonly) {
            auto <- mget(".trackAuto", ifnotfound=list(list(on=FALSE, last=-1)), envir=trackingEnv)[[1]]
            ## If this is an auto-update tracking env, get vars and files in sync,
            ## otherwise we presume the user knows what they're doing -- ignore
            ## out-of-sync things
            if (auto$on) {
                if (debug)
                    cat("Calling track.sync(envir = ", environmentName(envir), ")\n", sep="")
                if (stop.on.error)
                    track.sync(envir=envir, master="envir", full=TRUE, trackingEnv=trackingEnv)
                else
                    try(track.sync(envir=envir, master="envir", full=TRUE, trackingEnv=trackingEnv))
            }
            if (debug)
                cat("Calling track.flush(envir = ", environmentName(envir), ")\n", sep="")
            if (stop.on.error)
                track.flush(envir=envir, force=TRUE)
            else
                try(track.flush(envir=envir, force=TRUE))
        }
        if (keepVars) {
            for (var in tracked.vars) {
                objVal <- getTrackedVar(var, trackingEnv=trackingEnv, opt=opt)
                remove(list=var, envir=envir)
                assign(var, objVal, envir=envir)
            }
        } else {
            ## Can't remove the variables from a locked environment -- these
            ## should all be active bindings.  This situation can arise when
            ## a readonly tracking env is being detached, so we can hope that
            ## the active bindings are garbage collected when the tracked
            ## environment is removed from the search path.
            if (!environmentIsLocked(envir)) {
                remove(list=tracked.vars, envir=envir)
            }
            if (sessionEnd && identical(envir, globalenv())) {
                ## Would be nice to print this message here, but
                ## printing messages from .Last is tricky because of the order
                ## of how things appear to user: first comes the question, then
                ## if the answer to "save?" is not cancel, .Last is run, and
                ## save.image() called.
                ## cat("All tracked objects removed from the global environment (pos=1) -- saying 'y' to save will only save untracked objects\n")
            }
        }
        ## Set the tracking env pointer to NULL rather than removing it, in
        ## case 'envir' is locked.
        setTrackingEnv(envir, NULL, readonly=opt$readonly)
        if (exists(".trackAuto", envir=trackingEnv, inherits=FALSE)) {
            ## Should always be OK to remove a var from the trackingEnv,
            ## even when the db is readonly or locked -- it's the tracked
            ## env that gets locked, not the tracking env.
            if (FALSE && opt$readonly) {
                current <- get(".trackAuto", envir=trackingEnv, inherits=FALSE)
                current$on <- FALSE
                assign(".trackAuto", current, envir=trackingEnv)
            } else {
                remove(list=".trackAuto", envir=trackingEnv)
            }
        }
        ## Assign a marker variable so that a finalizer can see
        ## when this env is done with.
        assign(".trackingFinished", TRUE, envir=trackingEnv)
        ## If we created this env via track.attach() and it is now empty
        ## of everything except tracking reserved names, detach it.
        if (detach && exists(".trackingCreated", envir=envir, inherits=FALSE)
            && !identical(globalenv(), envir)) {
            vars <- ls(envir=envir, all.names=TRUE)
            pos <- match(environmentName(envir), search())
            vars <- vars[!isReservedName(vars)]
            if (length(vars)==0 || force.detach || environmentIsLocked(envir)) {
                if (!is.na(pos)) {
                    if (verbose)
                        cat("Removing", envname(envir), "from the search path\n")
                    detach(pos=pos)
                } else {
                    cat("Strange: can't find", envname(envir), "on the search path\n")
                }
            } else {
                warning("No longer tracking it, but can't remove ", envname(envir), " from the search path: still contains some variables: ", paste(vars[!isReservedName(vars)], collapse=", "), "\n", sep="")
            }
        }
    }
    return(invisible(NULL))
}
