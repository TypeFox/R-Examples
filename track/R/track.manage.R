track.remove <- function(expr, pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, all=FALSE, force=TRUE)
    trackedVarOp(if (!missing(expr)) substitute(expr), envir=envir, list=list, pattern=pattern, glob=glob, all=all, op="remove", who="track.remove()", force=force)

track.save <- function(expr, pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, all=missing(expr) && missing(list) && missing(pattern) && missing(glob))
    trackedVarOp(if (!missing(expr)) substitute(expr), envir=envir, list=list, pattern=pattern, glob=glob, all=all, op="save", resave=FALSE, who="track.save()")

track.resave <- function(expr, pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, all=missing(expr) && missing(list) && missing(pattern) && missing(glob))
    trackedVarOp(if (!missing(expr)) substitute(expr), envir=envir, list=list, pattern=pattern, glob=glob, all=all, op="save", resave=TRUE, who="track.resave()")

track.flush <- function(expr, pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, all=missing(expr) && missing(list) && missing(pattern) && missing(glob), force=FALSE)
    trackedVarOp(if (!missing(expr)) substitute(expr), envir=envir, list=list, pattern=pattern, glob=glob, all=all, op="flush", who="track.flush()", force=force)

track.forget <- function(expr, pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, all=FALSE)
    trackedVarOp(if (!missing(expr)) substitute(expr), envir=envir, list=list, pattern=pattern, glob=glob, all=all, op="forget", who="track.forget()")

untrack <- function(expr, pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, all=FALSE, keep.in.db=FALSE)
    trackedVarOp(if (!missing(expr)) substitute(expr), envir=envir, list=list, pattern=pattern, glob=glob, all=all, op=if (keep.in.db) "lift" else "untrack", who="untrack()")

## perform some operation on tracked variables
trackedVarOp <- function(qexpr, pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, all=FALSE,
                         op=c("save", "flush", "forget", "remove", "untrack", "lift", "orphan"),
                         resave=FALSE, force=FALSE, who="?") {
    op <- match.arg(op)
    trackingEnv <- getTrackingEnv(envir)
    opt <- track.options(trackingEnv=trackingEnv)
    dir <- getTrackingDir(trackingEnv)
    dataDir <- getDataDir(dir)
    fileMap <- getFileMapObj(trackingEnv)
    unsaved <- getUnsavedObj(trackingEnv, NULL)
    objSummary <- getObjSummary(trackingEnv, opt=opt)
    if (!is.null(qexpr)) {
        if (is.name(qexpr) || is.character(qexpr)) {
            objName <- as.character(qexpr)
        } else {
            stop("expr argument to ", who, " must be a quoted or unquoted variable name")
        }
        list <- c(objName, list)
    }
    if (!is.null(pattern) || !is.null(glob) || all) {
        if (!is.null(list))
            stop("cannot use expr= or list= at the same time as pattern=, glob=, or all=TRUE")
        list <- track.status(envir=envir, pattern=pattern, glob=glob, file.status=FALSE,
                             what=if (op=="save") "unsaved" else "tracked", tracked=TRUE, all.names=TRUE)
    }
    all.objs <- ls(envir=envir, all.names=TRUE)
    if (length(list)) {
        isTracked <- objIsTracked(list, envir, trackingEnv, all.objs)
        if (!force && !all(isTracked)) {
            cat("The following objects are not tracked: ",
                paste("'", list[which(!isTracked)[seq(len=min(3,sum(!isTracked)))]], "'", sep="", collapse=", "),
                if (sum(!isTracked) > 3) ", ...",
                "\n", sep="")
            list <- list[!isTracked]
        }
    }
    i <- which(isReservedName(list))
    if (length(i)) {
        warning("cannot ", op, " ", length(i), " vars (these are used in implementing tracking): ",
                paste("'", list[i], "'", sep="", collapse=", "))
        list <- list[-i]
    }
    quarantine.dir <- file.path(dataDir, "quarantine")
    needSaveFileMap <- FALSE
    needSaveObjSummary <- mget(".trackingSummaryChanged", envir=trackingEnv, ifnotfound=list(FALSE))[[1]]
    track.preremove.methods <- character(0)
    for (objName in list) {
        fileMapChanged <- FALSE
        objSummaryChanged <- FALSE
        # if (!exists(objName, envir=envir, inherits=FALSE)) { # don't use exists() because it ...
        if (!force && !is.element(objName, all.objs)) {
            warning("'", objName, "' does not exist in ", envname(envir))
            next
        }
        if (is.element(op, c("untrack", "lift")) && !force && !bindingIsActive(objName, envir)) {
            warning("cannot ", op, " tracked var '", objName, "' because it is not a properly tracked variable (not an active binding in ", envname(envir), ")")
            next
        }
        fileMapPos <- match(objName, names(fileMap))
        if (is.na(fileMapPos)) {
            if (objIsTracked(objName, envir, trackingEnv, all.objs)) {
                warning("cannot ", op, " tracked var '", objName, "' because it is not in the file map in ", envname(envir))
            } else {
                if (op=="remove")
                    remove(list=objName, envir=envir, inherits=FALSE)
            }
            next
        }
        file <- paste(fileMap[fileMapPos], opt$RDataSuffix, sep=".")
        filePath <- file.path(dataDir, file)
        ##
        ## First task is to modify envir if necessary
        ##
        if (op=="remove") {
            ## remove the variable from envir
            if (is.element(objName, all.objs)) {
                if (objName %in% rownames(objSummary)) {
                    ## See if we need to call a track.preremove() method.
                    ## Don't want to actually fetch the object to do this unless it is necessary.
                    cls <- objSummary[objName, "class"]
                    found <- FALSE
                    for (cl in strsplit(cls, ",")[[1]]) {
                        meth <- paste("track.preremove", cl, sep=".")
                        if (meth %in% track.preremove.methods) {
                            found <- TRUE
                            break
                        } else if (exists(meth, mode="function")) {
                            found <- TRUE
                            track.preremove.methods <- c(track.preremove.methods, meth)
                            break
                        }
                    }
                    if (found) {
                        objVal <- get(objName, envir=envir, inherits=FALSE)
                        res <- try(track.preremove(objVal, objName, envir), silent=TRUE)
                    }
                }
                remove(list=objName, envir=envir)
            }
        } else if (is.element(op, c("untrack", "lift"))) {
            ## Fetch the value and assign it in envir
            ## Don't run any load callback on the object
            if (exists(objName, envir=trackingEnv, inherits=FALSE)) {
                objVal <- get(objName, envir=trackingEnv, inherits=FALSE)
            } else {
                tmpenv <- new.env(parent=emptyenv())
                load.res <- load(filePath, envir=tmpenv)
                if (length(load.res)<1 || load.res[1] != objName) {
                    warning("file '", filePath, "' did not contain obj '", objName,
                            "' as its first object - moving this file to ", quarantine.dir)
                    if (!dir.exists(quarantine.dir))
                        dir.create(quarantine.dir)
                    file.rename(filePath, file.path(dataDir, "quarantine", file))
                } else {
                    if (length(load.res)>1) {
                        warning("ignoring other objects in file '", filePath, "'")
                        remove(list=load.res[-1], envir=tmpenv)
                    }
                    objVal <- get(objName, envir=tmpenv, inherits=FALSE)
                    remove(list=objName, envir=tmpenv)
                }
            }
            # need to remove the active binding from envir before saving the ordinary variable
            remove(list=objName, envir=envir)
            assign(objName, objVal, envir=envir)
        } else if (!is.element(op, c("save", "flush", "forget"))) {
            stop("what ", op, "???")
        }
        ##
        ## Second task is to remove files and modify trackingEnv if necessary
        ##
        if (is.element(op, c("remove", "untrack"))) {
            fileMap <- fileMap[-fileMapPos]
            fileMapChanged <- TRUE
            if (opt$debug >= 2)
                cat('trackedVarOp:', op, 'removing', paste(objName, collapse=', '), 'from trackingEnv\n')
            if (exists(objName, envir=trackingEnv, inherits=FALSE))
                remove(list=objName, envir=trackingEnv)
            if (file.exists(filePath)) {
                rm.res <- try(file.remove(filePath), silent=TRUE)
                if (is(rm.res, "try-error"))
                    warning("could not remove file for tracked var '", objName, "'")
            }
            i <- match(objName, rownames(objSummary))
            if (!is.na(i)) {
                objSummary <- objSummary[-i,,drop=FALSE]
                objSummaryChanged <- TRUE
            }
        } else if (is.element(op, c("save", "flush", "forget", "lift"))) {
            if (is.element(op, c("flush", "save", "lift")) && exists(objName, envir=trackingEnv, inherits=FALSE)
                && (resave || is.element(objName, unsaved))) {
                save.res <- try(save(list=objName, envir=trackingEnv, file=filePath,
                                     compress=opt$compress, compression_level=opt$compression_level), silent=TRUE)
                if (is(save.res, "try-error"))
                    stop("could not save '", objName, "' in ", filePath, ": fix file problem and try again")
            }
            if (is.element(op, c("flush", "forget", "lift")) && exists(objName, envir=trackingEnv, inherits=FALSE)
                && (force || !(op=="flush" && is.element(objName, opt$alwaysCache)))) {
                if (opt$debug >= 2)
                    cat('trackedVarOp:', op, 'removing', paste(objName, collapse=', '), 'from trackingEnv\n')
                remove(list=objName, envir=trackingEnv)
            }
        } else {
            stop("what ", op, "???")
        }
        ## every operation we do has the effect of unsetting the 'unsaved' bit
        if (!is.na(i <- match(objName, unsaved))) {
            unsaved <- unsaved[-i]
            assign(".trackingUnsaved", unsaved, envir=trackingEnv)
        }
        if (fileMapChanged) {
            needSaveFileMap <- TRUE
            assign(".trackingFileMap", fileMap, envir=trackingEnv)
        }
        if (objSummaryChanged) {
            needSaveObjSummary <- TRUE
            assign(".trackingSummary", objSummary, envir=trackingEnv)
        }
    }
    save.res <- NULL
    ## Note that we never write the "unsaved" list out to a file -- it just stays in memory
    if ((needSaveFileMap || resave) && !opt$readonly)
        writeFileMapFile(fileMap, trackingEnv, dataDir, FALSE)
    if ((needSaveObjSummary || resave) && !opt$readonly) {
        assign(".trackingSummaryChanged", TRUE, envir=trackingEnv)
        save.res <- saveObjSummary(trackingEnv, opt, dataDir)
        if (!is(save.res, "try-error"))
            assign(".trackingSummaryChanged", FALSE, envir=trackingEnv)
    }
    if (is(save.res, "try-error")) {
        warning("unable to save .trackingSummary in ", file.path(dataDir), ": if this message appears repeatedly, fix problem and run track.resave()")
    }
    return(invisible(list))
}

track.preremove <- function(obj, objName, envir, ...) UseMethod("track.preremove", obj)
