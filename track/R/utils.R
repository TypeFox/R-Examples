## Various internal utility functions for trackObjs

## Design:
##
## There is a trackingEnv variable in the tracked
## environment (e.g., .GlobalEnv)
##
## Tracked vars are active bindings, whose function refers to their
## name and trackingEnv
##
## The trackingEnv has several pieces of metadata:
##   * .trackingSummary
##   * .trackingFileMap (contains all objects)
##   * .trackingUnsaved (contains names of unsaved objects, only used if writeToDisk==FALSE)
##   * .trackOptions (list(writeToDisk=TRUE/FALSE, cache=TRUE/FALSE)) (default TRUE/TRUE)
##   * .trackingDir (absolute path -- don't depend on getwd() -- can be changed
## The first two are stored on disk, and written back out whenever
## they are changed.
## .trackingUnsaved is only stored in memory.
## When the database is initially read, .trackingFileMap and .trackingSummary
## might be created.  Note that object names are stored in RData files,
## so we can recreate these objects if they are missing, and check
## consistency when we read them.
##
## Functions: ( * = to be written)
##   track.start(dir=trackingDir)
##   track.dir()
##   track.stop()
##   track(var) track(var <- value) track(list=...) track(all=TRUE)
##   track.status(): return a list of track env, track dir, and tracked, untracked, and untrackable vars
##   untrack(var) untrack(list=) untrack(all=TRUE)
##   track.remove(var) track.remove(list=) track.remove(all=TRUE) (call this track.remove() ?)
##   track.summary()
##
##   tracked(): return all tracked variables
##   untracked(): return all untracked, but trackable variables
##   untrackable(): return all untrackable variables
##   track.unsaved(): return all tracked but unsaved variables
##   track.orphaned(): return variables that exist in the tracking dir, but have no binding
##   track.masked(): return variables that exist in the tracking dir, but have an ordinary binding
##
##   track.flush(var, list, all): write unsaved variables to disk, and remove from memory
##   track.save(var, list, all, book): write unsaved variables to disk
##   track.forget(var, list, all): delete cached versions without saving to file
##   track.rescan(): reload info from disk (forget all cached vars, might remove some no-longer existing tracked vars)
##   track.rebuild(dir=trackingDir)
##
## Don't need to use delayedAssign if use makeActiveBinding
## vars that have a binding, and no corresponding var in trackingEnv
## means the var must be out on disk.  However, can still use
## the same structure as a ddp from g.data().
##
## database of object data is stored in .trackingSummary
## contains the following columns:
##   name
##   class
##   mode
##   extent
##   length
##   size
##   time: modified
##   time: created
##   time: accessed
##   ES: sessions alive
##   SA: num accesses this session
##   SW: num writes this session
##   PA: total accesses before this session
##   PW: total writes before this session
##
## Also need to have mapping of object names -> file names:
## store this in .trackingFileMap.
## Any names that contain upper case letters or don't conform
## to starting with alphanum and then continuing with alphanum . _
## are given a filename like .XXXX where XXXX is a number.
##
## When reconstructing trackingSummary, use file.info() to get file
## creation, modification and last access times.
##
## Rules for names of RData files:
##    'simple' names satisfy the following conditions:
##       * contain lowercase letters, digits, '.', "_"
##       * start with a lowercase letter
##       * are 55 characters or less
##       * are not one of 'con', 'prn', 'aux', 'nul', 'com1' .. 'com9', 'lpt1' .. 'lpt9'
##    * when an object has a filename that does not begin with a '_',
##      the filename and the object name must be simple and identical
##    * objects that do not have a simple name must be stored in
##      files with names beginning with a '_'

env.is.tracked <- function(pos=1, envir=as.environment(pos)) {
    # .trackingEnv could be present and NULL for a non-tracked environment
    env <- mget(".trackingEnv", ifnotfound=list(NULL), envir=envir)[[1]]
    return(!is.null(env)
           && exists(".trackingPid", envir=env, inherits=FALSE)
           && identical(get(".trackingPid", envir=env, inherits=FALSE), Sys.getpid()))
}

getDataDir <- function(trackingDir) {
    ## If we wanted to store saved RData files in a subdirectory called 'data',
    ## we would use this here:
    ## return(file.path(dir, "data"))
    trackingDir
}

tracked.envs <- function(envirs=search()) {
    ## returns environment names, not environments!
    i <- sapply(envirs, function(envir) env.is.tracked(envir=as.environment(envir)))
    envirs[i]
}


# Create or update a row for the summary data frame
summaryRow <- function(name, opt, sumRow=NULL, obj=NULL, file=NULL, change=FALSE, times=NULL, accessed=TRUE) {
    if (opt$use.fake.Sys.time)
        tt <- fake.Sys.time()
    else
        tt <- Sys.time()
    new <- FALSE
    if (is.null(sumRow)) {
        new <- TRUE
        sumRow <- data.frame(row.names=name, class="", mode="", extent="",
                             length=length(obj), size=as.double(object.size(obj)),
                             modified=tt, created=tt, accessed=tt,
                             A=as.integer(1), ES=as.integer(1), SA=as.integer(0),
                             SW=as.integer(0), PA=as.integer(0), PW=as.integer(0),
                             cache=NA, stringsAsFactors=FALSE)
        if (!is.null(file) && file.exists(file)
            && !is(info <- try(file.info(file), silent=TRUE), "try-error")) {
            if (length(info$mtime)==1 && !is.na(info$mtime))
                sumRow$modified <- info$mtime
            if (length(info$ctime)==1 && !is.na(info$ctime))
                sumRow$created <- info$ctime
            if (length(info$atime)==1 && !is.na(info$atime))
                sumRow$accessed <- info$atime
        }
    }

    if (change || new) {
        cl <- class(obj)
        ## Neither the first or last element of the returned value of 'class' is
        ## clearly the most useful class to note, so include all classes as a comma
        ## separated string.
        ## E.g., class of an object returned by 'glm()': "glm" "lm"
        ## while class of an object returned by Sys.time(): "POSIXt" "POSIXct"
        sumRow$class <- if (length(cl)==0) "?" else paste(cl, collapse=",")
        sumRow$mode <- mode(obj)
        if (new || (change && (is.na(sumRow$cache) || ! (sumRow$cache %in% c("fixedyes", "fixedno"))))) {
            sumRow$cache <- factor(ifelse(   is.element(name, opt$alwaysCache)
                                          || any(is.element(cl, opt$alwaysCacheClass)),
                                          "yes", "no"),
                                   levels=c("fixedno", "no", "yes", "fixedyes"))
        }
        l <- try(length(obj), silent=TRUE)
        if (is(l, "try-error"))
            sumRow$length <- NA
        else
            sumRow$length <- as.integer(l)
        d <- try(dim(obj), silent=TRUE)
        if (is(d, "try-error"))
            sumRow$extent <- "(error)"
        else if (is.numeric(d))
            sumRow$extent <- paste("[", paste(d, collapse="x"), "]", sep="")
        else
            sumRow$extent <- paste("[", format(sumRow$length), "]", sep="")
        if (is.list(obj))
            sumRow$extent <- paste("[", sumRow$extent, "]", sep="")
        sumRow$size <- object.size(obj)
        if (accessed)
            sumRow$modified <- tt
        sumRow$SW <- sumRow$SW + as.integer(1)
    } else if (accessed) {
        sumRow$accessed <- tt
        sumRow$SA <- sumRow$SA + as.integer(1)
    }
    if (!is.null(times)) {
        sumRow$modified <- times$mtime
        sumRow$created <- times$ctime
        sumRow$accessed <- times$atime
    }

    return(sumRow)
}

isReservedName <- function(objName)
    ## ".trackingEnv" is a reserved name to allow for storing the
    ## tracking env as an object in the tracked environment instead
    ## of an attribute on the tracked environment.
    return((regexpr("[\n\r]", objName) > 0)
           | is.element(objName, c(".trackingEnv", ".trackingDir", ".trackingFileMap",
                                   ".trackingUnsaved", ".trackingSummary",
                                   ".trackingSummaryChanged", ".trackingOptions",
                                   ".trackingPid", ".trackingCreated", ".trackingCacheMark",
                                   ".trackAuto", ".trackingFinished", ".trackingModTime")))

objIsTracked <- function(objNames, envir, trackingEnv, all.objs=ls(envir=envir, all.names=TRUE)) {
    if (length(objNames)==0)
        return(logical(0))
    fileMap <- getFileMapObj(trackingEnv)
    return(sapply(objNames, function(objName) {
        ## an object is already tracked if the following 2 conditions are met:
        ##   - it exists as an activing binding in envir
        ##   - there is an entry in the fileMap in the trackingEnv
        ## Don't use exists() because it gets the object
        ## if (!exists(objName, envir=envir, inherits=FALSE))
        if (!is.element(objName, all.objs))
            return(FALSE)
        if (!bindingIsActive(objName, envir))
            return(FALSE)
        if (!is.element(objName, names(fileMap)))
            return(FALSE)
        return(TRUE)
    }))
}

getTrackingDir <- function(trackingEnv) {
    dir <- NULL
    ## avoid partial name matching for attributes...
    dir <- mget(".trackingDir", ifnotfound=list(NULL), envir=trackingEnv)[[1]]
    if (is.null(dir))
        stop("trackingEnv ", envname(trackingEnv), " has no '.trackingDir' variable")
    if (!is.character(dir) || length(dir)!=1)
        stop("variable '.trackingDir' in trackingEnv ", envname(trackingEnv),
             " must be a length one character vector")
    return(dir)
}

setTrackingEnv <- function(trackedEnv, trackingEnv, readonly=FALSE) {
    ## When trackingEnv=NULL, this function should remove the tracking env
    ## if it can, otherwise set it to NULL.
    if (is.null(trackingEnv) && !environmentIsLocked(trackedEnv)) {
        if (exists(".trackingEnv", envir=trackedEnv, inherits=FALSE))
            remove(list=".trackingEnv", envir=trackedEnv)
    } else {
        assign(".trackingEnv", trackingEnv, envir=trackedEnv)
    }
    invisible(NULL)
}

getTrackingEnv <- function(trackedEnv, stop.on.not.tracked = TRUE) {
    env <- mget(".trackingEnv", ifnotfound=list(NULL), envir=trackedEnv)[[1]]
    if (is.null(env))
        if (stop.on.not.tracked)
            stop("env ", envname(trackedEnv), " is not tracked (has no '.trackingEnv' variable)")
        else
            return(NULL)
    if (!is.environment(env))
        stop("variable '.trackingEnv' in env ", envname(trackedEnv),
             " is not an environment")
    return(env)
}

getFileMapObj <- function(trackingEnv) {
    fileMap <- mget(".trackingFileMap", envir=trackingEnv, ifnotfound=list(NULL))[[1]]
    if (is.null(fileMap) || !is.character(fileMap))
        stop("no usable .trackingFileMap object in tracking env ", envname(trackingEnv), " - recommend using track.rebuild()")
    return(fileMap)
}

writeFileMapFile <- function(fileMap, trackingEnv, dataDir, assignObj=TRUE) {
    if (assignObj && is(try(assign(".trackingFileMap", fileMap, envir=trackingEnv), silent=TRUE), "try-error"))
        warning("failed to assign '.trackingFileMap' in ", envname(trackingEnv))
    if (length(fileMap)) {
        i <- order(names(fileMap))
        fileData <- paste(fileMap[i], ":", names(fileMap)[i], sep="")
    } else {
        fileData <- character(0)
    }
    ## open in binary mode so that we use just "\n" as a separator
    open.res <- (con <- file(file.path(dataDir, "filemap.txt"), open="wb"))
    if (is(open.res, "try-error")) {
        warning("failed to open ", file.path(dataDir, "filemap.txt"), " for writing: if this message appears repeatedly try to fix problem, then do 'track.resave()'")
        return(FALSE)
    }
    on.exit(close(con))
    save.res <- try(writeLines(text=fileData, con=con, sep="\n"), silent=TRUE)
    if (is(save.res, "try-error")) {
        warning("failed to save filemap.txt: if this message appears repeatedly try to fix problem, then do 'track.resave()'")
        return(FALSE)
    }
    return(TRUE)
}

readFileMapFile <- function(trackingEnv, dataDir, assignObj) {
    ## open in binary mode so that we use just "\n" as a separator
    open.res <- (con <- file(file.path(dataDir, "filemap.txt"), open="rb"))
    if (is(open.res, "try-error"))
        stop("failed to open \"", file.path(dataDir, "filemap.txt"), "\" for reading: try using track.rebuild()")
    on.exit(close(con))
    fileData <- try(readLines(con=con, n=-1), silent=TRUE)
    if (is(fileData, "try-error"))
        stop("failed to read file map data from \"", file.path(dataDir, "filemap.txt"), "\": try using track.rebuild()")
    ## Remove Windows line termination
    fileData <- gsub("\r", "", fileData)
    i <- regexpr(":", fileData, fixed=TRUE)
    if (any(i < 1))
        stop("file map contains invalid data (need a ':' in each line): \"", file.path(dataDir, "filemap.txt"), "\": try using track.rebuild()")
    fileMap <- substring(fileData, 1, i-1)
    names(fileMap) <- substring(fileData, i+1)
    if (assignObj && is(try(assign(".trackingFileMap", fileMap, envir=trackingEnv), silent=TRUE), "try-error"))
        warning("failed to assign '.trackingFileMap' in ", envname(trackingEnv))
    return(fileMap)
}

getObjSummary <- function(trackingEnv, opt, stop.if.not.found=TRUE, fromWhere=paste('tracking env ', envname(trackingEnv))) {
    objSummary <- mget(".trackingSummary", envir=trackingEnv, ifnotfound=list(NULL))[[1]]
    problem <- NULL
    if (is.null(objSummary))
        problem <- paste("no .trackingSummary object ", fromWhere,
                         " - recommend using track.rebuild()", sep='')
    else if (!is.data.frame(objSummary))
        problem <- paste(".trackingSummary object found in ", fromWhere,
                         " but is not a data frame - recommend using track.rebuild()", sep='')
    if (!is.null(problem))
        if (stop.if.not.found)
            stop(problem)
        else
            return(structure(problem, class='try-error'))
    ## The 'cache' column was added in version 1.03.
    ## If we read a .trackingSummary without it, just add it.
    if (!is.element("cache", names(objSummary)))
        objSummary$cache <- factor(rep(NA, nrow(objSummary)), levels=c("fixedno", "no", "yes", "fixedyes"))
    i <- is.na(objSummary$cache)
    objSummaryChanged <- FALSE
    ## If there are any NA values in 'cache', set them based on name and class.
    if (any(i)) {
        j <- is.element(rownames(objSummary)[i], opt$alwaysCache)
        if (any(j))
            objSummary[which(i)[j], "cache"] <- "yes"
        i <- is.na(objSummary$cache)
        if (any(i) & length(opt$alwaysCacheClasses)) {
            j <- sapply(strsplit(objSummary$class[i], ","), is.element, opt$alwaysCacheClasses)
            if (any(j))
                objSummary[which(i)[j], "cache"] <- "yes"
        }
        i <- is.na(objSummary$cache)
        if (any(i))
            objSummary[i, "cache"] <- "no"
        objSummaryChanged <- TRUE
    }
    if (objSummaryChanged) {
        assign(".trackingSummary", objSummary, envir=trackingEnv)
        assign(".trackingSummaryChanged", TRUE, envir=trackingEnv)
    }
    return(objSummary)
}

getUnsavedObj <- function(trackingEnv, notfound=character(0))
    mget(".trackingUnsaved", envir=trackingEnv, ifnotfound=list(notfound))[[1]]

envname <- function(envir) {
    # Produce a 1-line name for the environment.
    # Use the name it has on the search() list, if possible.
    # This is simpler now that R has the environmentName() function
    n <- environmentName(envir)
    if (n!="") {
        return(paste("<env ", n, ">", sep=""))
    }
    return(capture.output(print(envir))[1])
}

notyetdone <- function(msg) cat("Not yet done: ", msg, "\n", sep="")

isSimpleName <- function(objName) {
    return(nchar(objName)<=55
           && regexpr("^[[:lower:]][._[:digit:][:lower:]]*$", objName, perl=TRUE)==1
           && regexpr("^(prn|aux|con|nul|com[1-9]|lpt[1-9])(\\.|$)", objName)<0)
}

makeObjFileName <- function(objName, fileNames) {
    ## check if we can use objName as a filename
    if (isSimpleName(objName))
        return(objName)
    ## fileNames is a list of vectors of file names that are already used.
    ## Generate a filename of the form _NNN (NNN is a number without a leading zero)
    ## work out what numbers have been used
    if (is.list(fileNames))
        used <- unlist(use.names=FALSE, lapply(fileNames, function(x) as.numeric(substring(grep("^_[1-9][0-9]*$", x, value=TRUE), 2))))
    else
        used <- as.numeric(substring(grep("^_[1-9][0-9]*$", fileNames, value=TRUE), 2))
    used <- unique(used)
    used <- sort(c(0, used, length(used)+2))
    ## find the first gap
    i <- used[c(diff(used)!=1, TRUE)][1] + 1
    if (is.element(i, used))
        stop("algorithm failure!")
    file <- paste("_", i, sep="")
    return(file)
}

setTrackedVar <- function(objName, value, trackingEnv, opt=track.options(trackingEnv=trackingEnv), times=NULL, file=NULL, doAssign=TRUE) {
    ## 'file' should be without path or suffix
    if (opt$readonly)
        stop("variable '", objName, "' cannot be changed -- it is in a readonly tracking environment")
    ## Set the tracked var, and write it to disk if required
    if (opt$debug)
        cat("setting tracked var '", objName, "' in ", envname(trackingEnv), "\n", sep="")
    ## Need to assign it in the tracking env, because save() requires
    ## an object in an env. Maybe we could skip this step when cache=FALSE,
    ## but then we'd need to try to find it in its original location (if any).
    ## No longer true: could remove doAssign arg -- in one use case, the var
    ## is already in the env, so don't need the assign.
    ##
    ## Robustness: what to do if the assign fails?
    if (doAssign)
        assign(objName, value, envir=trackingEnv)
    ## Find the directory where we are saving, and create subdirs if necessary
    dir <- getTrackingDir(trackingEnv)
    for (d in unique(c(dir, getDataDir(dir))))
        if (!dir.exists(d))
            dir.create(d)
    ## Work out the name of the file to use for this var
    if (is.null(file)) {
        fileMap <- getFileMapObj(trackingEnv)
        fileMapChanged <- FALSE
        isNew <- FALSE
        if (is.na(file <- fileMap[match(objName, names(fileMap))])) {
            file <- makeObjFileName(objName, fileMap)
            fileMap[objName] <- file
            fileMapChanged <- TRUE
            isNew <- TRUE
        }
        if (fileMapChanged) {
            ##  always write a changed file map back out to disk
            writeFileMapFile(fileMap, trackingEnv=trackingEnv, dataDir=getDataDir(dir))
        }
    }
    fullFile <- NULL
    if (opt$writeToDisk && !is.element("eotPurge", opt$cachePolicy)) {
        ##  the value of 'file' is the base of the filename -- work out the full pathname
        fullFile <- file.path(getDataDir(dir), paste(file, opt$RDataSuffix, sep="."))
        if (opt$debug)
            cat("saving '", objName, "' to file ", fullFile, "\n", sep="")
        save.res <- try(save(list=objName, file=fullFile, envir=trackingEnv,
                             compress=opt$compress, compression_level=opt$compression_level), silent=TRUE)
        if (!is(save.res, "try-error")) {
            if (opt$debug >= 2)
                cat('setTrackedVar: removing', paste(objName, collapse=', '), 'from trackingEnv\n')
            if (!opt$cache && !is.element(objName, opt$alwaysCache))
                remove(list=objName, envir=trackingEnv)
            unsaved <- getUnsavedObj(trackingEnv)
            if (length(unsaved) && is.element(objName, unsaved))
                if (length(unsaved)>1)
                    assign(".trackingUnsaved", setdiff(unsaved, objName), envir=trackingEnv)
                else
                    remove(".trackingUnsaved", envir=trackingEnv)
        } else {
            stop("failed to save obj '", objName, "' in file '", fullFile,
                 "' (", as.character(save.res), ") - '", objName,
                 "' is currently in env ", envname(trackingEnv), " but is not saved to disk",
                 " (suggestion: fix disk problems then do track.save())")
        }
    } else {
        ## mark the object as unsaved
        unsaved <- getUnsavedObj(trackingEnv)
        if (!is.element(objName, unsaved))
            assign(".trackingUnsaved", sort(c(objName, unsaved)), envir=trackingEnv)
    }
    if (opt$maintainSummary) {
        objSummary <- getObjSummary(trackingEnv, opt=opt)
        if (!is.data.frame(objSummary)) {
            warning(".trackingSummary in ", envname(trackingEnv), " is not a data.frame: not updating summary; run track.rebuild()")
        } else {
            if (is.element(objName, rownames(objSummary))) {
                sumRow <- summaryRow(objName, opt=opt, sumRow=objSummary[objName, , drop=FALSE], obj=value,
                                        file=NULL, change=TRUE, times=times)
            } else {
                ## Don't have a row in the summary for this object.
                if (!isNew) {
                    if (is.null(fullFile))
                        fullFile <- file.path(getDataDir(dir), paste(file, opt$RDataSuffix, sep="."))
                    if (!file.exists(fullFile)) {
                        ## This can happen when an object is created via track.ff() -- want to avoid noise
                        ## so don't issue this warning.
                        ## warning("file '", fullFile, "' does not yet exist (for obj '", objName, "') - will create it when needed")
                        fullFile <- NULL
                    }
                }
                ## summaryRow will get times from the file if we supply it -- only
                ## want to do this in the exceptional case that this object was not
                ## new, but didn't have a summary row.
                sumRow <- summaryRow(objName, opt=opt, obj=value, file=if (!isNew) fullFile else NULL,
                                     change=TRUE, times=times)
                ## If this is not a new object, record that the times are not accurrate.
                if (!isNew)
                    sumRow$A <- 0
            }
            objSummary[objName, ] <- sumRow
            assign.res <- try(assign(".trackingSummary", objSummary, envir=trackingEnv), silent=TRUE)
            if (is(assign.res, "try-error")) {
                warning("unable to assign .trackingSummary back to tracking env on ",
                        envname(trackingEnv), ": ", assign.res)
            } else {
                assign(".trackingSummaryChanged", TRUE, envir=trackingEnv)
                if (opt$writeToDisk && !is.element("eotPurge", opt$cachePolicy)) {
                    save.res <- saveObjSummary(trackingEnv, opt=opt, dataDir=getDataDir(dir))
                    if (is(save.res, "try-error"))
                        warning("unable to save .trackingSummary to ", dir)
                    else
                        assign(".trackingSummaryChanged", FALSE, envir=trackingEnv)
                }
            }
        }
    }
    ## return the value of the variable
    return(invisible(value))
}

getTrackedVar <- function(objName, trackingEnv, opt=track.options(trackingEnv=trackingEnv)) {
    ## Get the tracked var if it exists, or load it from disk if it doesn't
    ## Special case when getting variable '.Last' -- don't want to stop with
    ## an error because that can cause an infinite loop, so if can't retrieve
    ## that var, then print a warning and return NULL
    if (!is.character(objName) || length(objName)!=1)
        stop("objName must be a length one character vector")
    if (exists(objName, envir=trackingEnv, inherits=FALSE)) {
        if (opt$debug)
            cat("getTrackedVar: '", objName, "' from cached obj in ", envname(trackingEnv), "\n", sep="")
        value <- get(objName, envir=trackingEnv, inherits=FALSE)
        dir <- NULL
        fullFile <- NULL
    } else {
        if (opt$debug)
            cat("getting tracked var '", objName, "' from file\n", sep="")
        fileMap <- getFileMapObj(trackingEnv)
        file <- fileMap[match(objName, names(fileMap))]
        if (is.na(file))
            if (objName=='.Last') {
                warning("'", objName, "' does not exist in tracking env ", envname(trackingEnv),
                     " (has no entry in .trackingFileMap); returning NULL")
                return(NULL)
            } else {
                stop("'", objName, "' does not exist in tracking env ", envname(trackingEnv),
                     " (has no entry in .trackingFileMap)")
            }
        dir <- getTrackingDir(trackingEnv)
        fullFile <- file.path(getDataDir(dir), paste(file, opt$RDataSuffix, sep="."))
        if (!file.exists(fullFile)) {
            if (objName=='.Last') {
                warning("file '", fullFile, "' does not exist (for obj '", objName, "'); returning NULL")
                return(NULL)
            } else {
                stop("file '", fullFile, "' does not exist (for obj '", objName, "')")
            }
        }
        tmpenv <- new.env(parent=emptyenv())
        ## hopefully, no actual copies are made while loading the object into tmpenv,
        ## then copying it to envir and returning it as the value of this function
        load.res <- load(fullFile, envir=tmpenv)
        ## call the load callback
        if (length(load.res)<1 || load.res[1] != objName)
            stop("file '", fullFile, "' did not contain obj '", objName, "' as its first object")
        if (length(load.res)>1)
            warning("ignoring ", length(load.res)-1, " other objects in file '", fullFile, "'")
        value <- get(objName, envir=tmpenv, inherits=FALSE)
        if (opt$cache)
            assign(objName, value, envir=trackingEnv)
        ## Call the load hooks based on class
        if (is.element('ff', class(value))) {
            fff <- attr(attr(value, 'physical'), 'filename')
            ffd <- dirname(fff)
            ffd2 <- dirname(ffd)
            fff1 <- basename(fff)
            fff2 <- basename(ffd)
            ## If the ff filename looks like .../ff/<varfile>.ff, and
            ## the file <trackingDir>/ff/<varfile>.ff exists, then
            ## switch to using that file.
            if (fff2 == 'ff' && fff1 == paste(file, '.ff', sep='')) {
                fffr <- file.path(dir, 'ff', fff1)
                if (file.exists(fffr) && fffr != fff) {
                    attr(attr(value, 'physical'), 'orig.filename') <- fff
                    attr(attr(value, 'physical'), 'filename') <- fffr
                }
            }
        }
    }
    ## Update the object table (object characteristics, accesses)
    if (opt$maintainSummary && opt$recordAccesses) {
        objSummary <- getObjSummary(trackingEnv, opt=opt)
        if (!is.data.frame(objSummary)) {
            warning(".trackingSummary in ", envname(trackingEnv), " is not a data.frame: not updating objSummary; run track.rebuild()")
        } else {
            if (is.element(objName, rownames(objSummary))) {
                sumRow <- summaryRow(objName, opt=opt, sumRow=objSummary[objName, , drop=FALSE], obj=value,
                                        file=NULL, change=FALSE, times=NULL)
            } else {
                if (is.null(fullFile)) {
                    fileMap <- getFileMapObj(trackingEnv)
                    file <- fileMap[match(objName, names(fileMap))]
                    if (is.na(file)) {
                        warning("'", objName, "' does not exist in tracking env ", envname(trackingEnv),
                             " (has no entry in .trackingFileMap)")
                    } else {
                        dir <- getTrackingDir(trackingEnv)
                        fullFile <- file.path(getDataDir(dir), paste(file, opt$RDataSuffix, sep="."))
                        if (!file.exists(fullFile)) {
                            warning("file '", fullFile, "' does not exist (for obj '", objName, "')")
                            fullFile <- NULL
                        }
                    }
                }
                sumRow <- summaryRow(objName, opt=opt, obj=value, file=fullFile, change=FALSE, times=NULL)
                ## Since this is not a new object, but this we are creating a new summary
                ## row for it, record that the times are not accurrate
                sumRow$A <- 0
            }
            objSummary[objName, ] <- sumRow
            assign.res <- try(assign(".trackingSummary", objSummary, envir=trackingEnv))
            if (is(assign.res, "try-error")) {
                warning("unable to assign .trackingSummary back to tracking env on ", envname(trackingEnv))
            } else {
                ## only makes sense to save() .trackingSummary if we were able to assign it
                assign(".trackingSummaryChanged", TRUE, envir=trackingEnv)
                if (opt$alwaysSaveSummary && !is.element("eotPurge", opt$cachePolicy)) {
                    if (is.null(dir))
                        dir <- getTrackingDir(trackingEnv)
                    save.res <- saveObjSummary(trackingEnv, opt=opt, dataDir=getDataDir(dir))
                    if (is(save.res, "try-error"))
                        warning("unable to save .trackingSummary to ", dir, ": ", save.res)
                    else
                        assign(".trackingSummaryChanged", FALSE, envir=trackingEnv)
                }
            }
        }
    }
    return(value)
}

quietNormalizePath <- function(path, winslash='/', mustWork=FALSE) {
    if (isTRUE(mustWork))
        return(normalizePath(path, winslash='/', mustWork=mustWork))
    newPath <- try(normalizePath(path, winslash='/', mustWork=mustWork), silent=TRUE)
    if (inherits(newPath, 'try-error'))
        newPath <- path
    allComp <- strsplit(gsub('[/\\\\]', winslash, newPath), winslash)[[1]]
    newComp <- allComp
    j <- 2
    for (i in seq(2, len=length(allComp)-1)) {
        if (allComp[i]=='.') {
            newComp <- newComp[-j]
        } else if (allComp[i]=='..') {
            newComp <- newComp[-c(j, j-1)]
            j <- max(1, j-1)
        } else {
            j <- j+1
        }
    }
    paste0(allComp, collapse=winslash)
}

find.relative.path <- function(path, file) {
    ## Find a way to express file as a relative path to path
    path <- quietNormalizePath(path, winslash='/')
    file <- quietNormalizePath(file, winslash='/')
    if (.Platform$OS.type == "windows") {
        path <- gsub("\\", "/", path, fixed=TRUE)
        file <- gsub("\\", "/", file, fixed=TRUE)
    }
    path.comp <- strsplit(path, split="/", fixed=TRUE)[[1]]
    file.comp <- strsplit(file, split="/", fixed=TRUE)[[1]]
    ## Windows doesn't normalize drive letters to all caps,
    ## which can result in a mismatch
    if (.Platform$OS.type == "windows") {
        if (length(grep("^[a-zA-Z]:$", path.comp[1])))
            path.comp[1] <- casefold(path.comp[1], upper=FALSE)
        if (length(grep("^[a-zA-Z]:$", file.comp[1])))
            file.comp[1] <- casefold(file.comp[1], upper=FALSE)
    }
    i <- 1
    while (i <= min(length(path.comp), length(file.comp))
           && path.comp[i] == file.comp[i])
        i <- i+1
    if (i>1) {
        file.rel <- file.comp[seq(to=length(file.comp), len=length(file.comp)-(i-1))]
        ## path.use <- path.comp[seq(from=1, len=i-1)]
        if (i <= length(path.comp))
            file.rel <- c(rep("..", length(path.comp)-(i-1)), file.rel)
        if (length(file.rel)==0)
            file.rel <- "."
        return(paste(file.rel, collapse=.Platform$file.sep))
    } else {
        return(paste(file.comp, collapse=.Platform$file.sep))
    }
}

exclude.from.tracking <- function(objName, objClasses=NULL, opt) {
    exclude <- rep(FALSE, length(objName))
    if (length(objName) > 1 && length(objClasses) && (!is.list(objClasses) || length(objClasses)!=length(objName)))
        stop("objClasses must be list the same length as objName when objName is a vector")
    if (length(objName)==1 && !is.list(objClasses) && length(objClasses))
        objClasses <- list(objClasses)
    if (length(opt$autoTrackExcludeClass) && length(objClasses))
        exclude <- sapply(objClasses, is.element, opt$autoTrackExcludeClass)
    for (re in opt$autoTrackExcludePattern)
        exclude <- exclude | regexpr(re, objName) >= 1
    exclude
}

dir.exists <- function(dir) {
    ## Need this because sometimes directories on network drives
    ## under windows aren't seen by file.exists().  This can be
    ## true even if they have files in them.  In such cases the
    ## tests for existence for the files in them seem to work fine.
    ## However, file.access() seems more reliable, so try both.
    file.exists(dir) || (file.access(dir, mode=0)==0)
}

## Use these to override Sys.time() for testing purposes (to get Sys.time() to return
## reproducible, steadily incrementing "times".  This  just counts 1 second forward
## from a fixed starting time each time it is called.
## Only use these for testing because they introduce overhead for Sys.time().

fake.Sys.time <- function() {
    fake.Sys.time.counter <- get("fake.Sys.time.counter", envir=fake.Sys.time.env)
    assign("fake.Sys.time.counter", fake.Sys.time.counter + 1, envir=fake.Sys.time.env)
    return(fake.Sys.time.counter)
}

fake.Sys.time.env <- as.environment(list(fake.Sys.time.counter=as.POSIXct("2001/01/01 09:00:00", tz="GMT")+1))

set.fake.Sys.time <- function(offset=1) {
    assign("fake.Sys.time.counter", as.POSIXct("2001/01/01 09:00:00", tz="GMT")+offset, envir=fake.Sys.time.env)
    return(invisible(get("fake.Sys.time.counter", envir=fake.Sys.time.env)))
}

saveObjSummary <- function(trackingEnv,
                           opt=track.options(trackingEnv=trackingEnv),
                           dataDir=getDataDir(getTrackingDir(trackingEnv)),
                           envir=trackingEnv) {
    ## Save the object summary out to the file system
    ## envir is the actual environment where the .trackingSummary object lives; it is
    ## usually the same as trackingEnv
    file <- file.path(dataDir, paste(".trackingSummary", opt$RDataSuffix, sep="."))
    if (!exists(".trackingSummary", envir=envir, inherits=FALSE))
        return(structure('saveObjSummary: .trackingSummary does not exist', class='try-error'))
    objSummary <- get(".trackingSummary", envir=envir, inherits=FALSE)
    ## pad just serves to make the saved tracking summary have different size on disk
    ## which helps with noticing when a tracking summary has changed, because if
    ## it changes within the same second, it will still have the same mod time.
    ## However, most of the time, we're not concerned with second-level accuracy.
    if (opt$debug) {
        pad <- attr(objSummary, 'pad')
        if (is.null(pad))
            pad <- raw(1)
        length(pad) <- (length(pad) + 17) %% 101
        attr(objSummary, 'pad') <- pad
        assign(".trackingSummary", objSummary, envir=envir)
    }
    save.res <- try(save(list=".trackingSummary", file=file, envir=envir, compress=FALSE), silent=TRUE)
    if (is(save.res, 'try-error'))
        attr(save.res, 'file') <- file
    if (identical(envir, trackingEnv)) {
        modTime <- file.info(file)
        if (!is.na(modTime$mtime))
            try(assign('.trackingModTime', modTime[,c('size','mtime')], envir=trackingEnv), silent=TRUE)
    }
    save.res
}

loadObjSummary <- function(trackingEnv,
                           opt=track.options(trackingEnv=trackingEnv),
                           dataDir=getDataDir(getTrackingDir(trackingEnv)),
                           stop.on.not.exists=FALSE) {
    # Load the object summary from the file system
    file <- file.path(dataDir, paste(".trackingSummary", opt$RDataSuffix, sep="."))
    tmpenv <- new.env(parent=emptyenv())
    if (!file.exists(file))
        if (stop.on.not.exists)
            stop("file storing object summary doesn't exist: ", file)
        else
            return(NULL)
    load.res <- try(load(file, envir=tmpenv), silent=TRUE)
    if (is(load.res, "try-error"))
        stop(file, " cannot be loaded -- for recovery see ?track.rebuild (",
             as.character(load.res), ")")
    if (length(load.res)!=1 || load.res != ".trackingSummary")
        stop(file, " does not contain just '.trackingSummary' -- for recovery see ?track.rebuild")
    modTime <- file.info(file)
    if (!is.na(modTime$mtime))
        try(assign('.trackingModTime', modTime[,c('size','mtime')], envir=trackingEnv), silent=TRUE)
    ## .trackingSummary has to exist in tmpenv because we just loaded it
    getObjSummary(tmpenv, opt=opt, fromWhere=file)
}

# used to store data around failures in an attempt to avoid infinite loops
track.private.env <- as.environment(list(standard.Last=FALSE, last.failed.get=character(100), last.failed.time=rep(Sys.time(), 100), last.failed.i=1))
