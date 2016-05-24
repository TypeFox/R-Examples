track.rebuild <- function(pos=1, envir=as.environment(pos), dir=NULL, fix=FALSE, level=c("missing", "all"), trust=c("unspecified", "environment", "db"), verbose=1, RDataSuffix=NULL, dryRun=TRUE, replace.wd=TRUE, use.file.times=TRUE) {
    ## Rebuild the fileMap and/or the trackingSummary in the tracking dir.
    ## Those objects generally contain redundant information, but they can be
    ## expensive to rebuild, because rebuilding them requires reading every
    ## saved file.  Furthermore, creation/modification/read times saved in
    ## the summary may be different to those on files.
    ##
    ## The strategy for rebuilding doesn't really depend on whether or not
    ## a tracking dir is active -- it's slightly easier if it's not,
    ## because then we don't have to worry about checking or installing
    ## active bindings.  However, in either case, we need to use what we
    ## can of an existing summary, replacing those parts that are not present.
    ##
    ## If fix=FALSE, the only files that can be changed are
    ## .trackingFileMap.rda and .trackingSummary.rda.  If there are other
    ## problems than missing, out-of-date, or unusable versions of these
    ## objects, this function will stop with an error without changing
    ## anything in the tracking database.
    ##
    ## If fix=TRUE, the following problems are fixed:
    ##     * unreadable .rda files are moved to the 'quarantine' directory
    ##     * unusable .rda files (multiple objects, illegal or duplicated names)
    ##       are moved to the 'quarantine' directory
    ##
    ## If dryRun=TRUE, nothing should be changed (but if fix=TRUE, fixes
    ## are reported but not made).
    ##
    ## This function is quite long and complex because it has to deal with
    ## multiple sources of potentially partial information.  Information can
    ## come from four sources:
    ##   (1) metadata in the database on disk (the filemap and summary files)
    ##   (2) metadata in the tracking environment (again, the filemap and summary files)
    ##   (3) info in the saved objects (object names, characteristics and modification times)
    ##   (4) info in the tracking environment (object names and characteristics)
    ## This function tries to preserve as much information as it can because
    ## it can be expensive to rebuild all the information for databases with
    ## hundreds or thousands of objects.
    ## Tricking situations:
    ##   * objects in the tracking environment that do not have corresponding file
    ##     or entry in the filemap or summary
    formatMsg <- function(msg) {
        # Reformat an error message for better appearance inline
        msg <- gsub("[ \t]*\n[ \t]*", " ", msg)
        msg <- gsub("[ \t]*$", "", msg)
        msg
    }

    chooseBestSummaryRow <- function(objName) {
        ## find the most recent row
        x <- envSummary[match(objName, row.names(envSummary), nomatch=0),,drop=FALSE]
        y <- fileSummary[match(objName, row.names(fileSummary), nomatch=0),,drop=FALSE]
        if (nrow(x) && nrow(y)) {
            y.newer <- (max(unlist(x[,c("modified", "created", "accessed")]), na.rm=TRUE)
                        < max(unlist(y[,c("modified", "created", "accessed")]), na.rm=TRUE))
            if (!is.na(y.newer) && y.newer)
                x <- y
        } else if (nrow(y)) {
            x <- y
        }
        if (nrow(x)==0)
            return(NULL)
        else
            return(x)
    }

    setNamedElt <- function(x, name, elt) {
        ## fill an empty space in x (an empty space has an empty name)
        ## or extend x if there are no empty spaces
        if (length(name)!=1 && length(name) != length(elt))
            stop("length(name) != length(elt)")
        for (j in seq(len=length(name))) {
            i <- match(name[j], names(x))
            if (is.na(i)) {
                if (is.null(names(x))) {
                    i <- 1
                } else {
                    i <- match("", names(x))
                    if (is.na(i))
                        i <- which(is.na(names(x)))[1]
                    if (is.na(i))
                        i <- length(x)+1
                }
            }
            if (length(name)==1)
                x[[i]] <- elt
            else
                x[[i]] <- elt[j]
            names(x)[i] <- name[j]
        }
        x
    }
    # replace backslash with forward slash because getAbsolutePath() always returns forward slashes
    wd <- gsub("\\\\", "/", getwd())
    abbrevWD <- function(path) {
        # abbrevWD replaces wd (working directory) in the path by "."
        # This makes the diagnostic output more concise and portable
        # across machines for matching in tests
        if (replace.wd)
            sub(wd, ".", path)
        else
            path
    }

    if (dryRun)
        cat("dryRun=TRUE: no changes will be actually be made\n")
    trackingEnv <- NULL
    if (!is.null(dir)) {
        dir <- getAbsolutePath(dir)
        if (!missing(envir) || !missing(pos))
            stop("must supply either 'dir' or 'envir' or 'pos'")
        ## Work out if this tracking dir is active
        ## Need to see if this directory is the one used on any
        ## tracking environments.  Get around the issue of aliasing
        ## (different names for the same directory) by creating a
        ## temporary file in the directory, and then looking for it.
        activeTracking <- FALSE
        if (!dir.exists(dir))
            stop("'", dir, "' does not exist")
        tmpfile <- tempfile(".rebuildTest", dir)
        tmpfilebase <- basename(tmpfile)
        if (!file.create(tmpfile))
            stop("unable to create temporary test file '", tmpfilebase, "' in dir '", dir, "'")
        on.exit(unlink(tmpfile))
        for (e.name in tracked.envs()) {
            e <- as.environment(e.name)
            e.dir <- getTrackingDir(getTrackingEnv(e))
            if (file.exists(file.path(e.dir, tmpfilebase))) {
                ## found the directory!
                envir <- e
                activeTracking <- TRUE
                trackingEnv <- getTrackingEnv(envir)
            }
        }
        unlink(tmpfile)
        on.exit()
        if (!activeTracking)
            envir <- NULL
    } else {
        if (!env.is.tracked(envir=envir))
            stop("env ", envname(envir), " is not tracked")
        activeTracking <- TRUE
        trackingEnv <- getTrackingEnv(envir)
        dir <- getTrackingDir(trackingEnv)
    }
    if (verbose)
        if (activeTracking)
            cat("Rebuilding an active tracking database\n")
        else
            cat("Rebuilding an inactive tracking database\n")
    dataDir <- getDataDir(dir)
    ## After here the following variables can be used:
    ##   'dir', 'dataDir'
    ##   'activeTracking'
    ##   'envir' (will be NULL if activeTracking==FALSE)
    ##   'trackingEnv'
    level <- match.arg(level)
    trust <- match.arg(trust)
    quarantineDir <- file.path(dataDir, "quarantine")
    ##
    ## Work out the options
    ## The trickiest option to work out is the RData file suffix, because
    ## we have to know it to find a saved .trackingOptions object, but this
    ## object may not be there.  So, first look for .trackingOptions.*, then
    ## then .trackingSummary.*, then data files.
    gopt <- getOption("global.track.options")
    if (length(gopt$RDataSuffixes)==0)
        gopt$RDataSuffixes <- c("rda", "RData")
    if (!is.character(gopt$RDataSuffixes))
        stop('getOption("global.track.options")$RDataSuffixes must be character data')
    if (any(!regexpr("^[[:alnum:]]+$", gopt$RDataSuffixes)))
        stop('getOption("global.track.options")$RDataSuffixes must consist of alpha-numeric characters only')
    if (!is.null(RDataSuffix)) {
        if (!is.character(RDataSuffix) || length(RDataSuffix)!=1)
            stop('RDataSuffix must be character data of length 1')
        if (any(!regexpr("^[[:alnum:]]+$", RDataSuffix)))
            stop('RDataSuffix must consist of alpha-numeric characters only')
    }
    if (length(gopt$RDataSuffixes)==1)
        suffixRegExp <- gopt$RDataSuffixes
    else
        suffixRegExp <- paste("(", paste(gopt$RDataSuffixes, collapse="|", sep=""), ")", sep="")
    if (!dir.exists(file.path(dataDir)))
        stop("dataDir does not exist: \"", dataDir, "\"")

    ## Try to work out the suffix being used
    ## First look for .trackingOptions file
    suffix <- NULL
    x <- list.files(path=dataDir, pattern=paste("^\\.trackingOptions\\.", suffixRegExp, "$", sep=""), all.files=TRUE)
    if (length(x)>1)
        stop("have multiple options files in '", dataDir, "': ", paste(x, collapse=", "))
    if (length(x)==1) {
        suffix <- sub(".*\\.", "", x)
        cat("Guessing RDataSuffix to be '", suffix, "' from name '", x, "'\n", sep="")
    }
    if (is.null(suffix)) {
        ## next look for .trackingSummary file
        x <- list.files(path=dataDir, pattern=paste("^\\.trackingSummary\\.", suffixRegExp, "$", sep=""), all.files=TRUE)
        if (length(x)>1)
            stop("have multiple summary files in '", dataDir, "': ", paste(x, collapse=", "))
        if (length(x)==1) {
            suffix <- sub(".*\\.", "", x)
            cat("Guessing RDataSuffix to be '", suffix, "' from name '", x, "'\n", sep="")
        }
    }
    if (is.null(suffix)) {
        ## next look for any files with possible RData suffix
        x <- list.files(path=dataDir, pattern=paste("^.*\\.", suffixRegExp, "$", sep=""), all.files=TRUE)
        if (length(x)>0) {
            suffix <- unique(sub(".*\\.", "", x))
            if (length(suffix)>1)
                stop("have files with multiple RData suffixes in '", dataDir, "': ", paste(x, collapse=", "))
            cat("Guessing RDataSuffix to be '", suffix, "' from names of various RData files: ",
                paste("'", x[max(length(x), 2)], "'", collapse=", ", sep=""),
                if (length(x)>2) "...", "\n", sep="")
        } else {
            if (!is.null(RDataSuffix)) {
                suffix <- RDataSuffix
                cat("Using RDataSuffix '", suffix, "' from supplied argument\n", sep="")
            } else {
                suffix <- gopt$RDataSuffixes[1]
                cat("Using RDataSuffix '", suffix, "' from options('global.track.options')\n", sep="")
            }
        }
    }
    if (!is.element(suffix, gopt$RDataSuffixes))
        stop("internal error: ended up with an illegal suffix?? (", suffix, ")")
    if (!is.null(RDataSuffix) && RDataSuffix != suffix)
        stop("suffix in use '", suffix, "' differs from supplied RDataSuffix ('", RDataSuffix, "')")

    if (activeTracking) {
        opt <- track.options(trackingEnv=trackingEnv)
        if (opt$readonly)
            stop("cannot rebuild a readonly tracking environment")
    } else {
        opt <- list(RDataSuffix=suffix)
        opt <- list()
        ## Read the options, using the suffix we found above
        file <- file.path(dataDir, paste(".trackingOptions", suffix, sep="."))
        if (file.exists(file)) {
            tmpenv <- new.env(parent=emptyenv())
            load.res <- try(load(file=file, envir=tmpenv), silent=TRUE)
            if (is(load.res, "try-error") || length(load.res)!=1 || load.res!=".trackingOptions") {
                warning(file, " does not contain a .trackingOptions object -- ignoring it and using system defaults")
            } else {
                opt <- get(".trackingOptions", envir=tmpenv, inherits=FALSE)
            }
        }
        opt <- track.options(old.options=opt, envir=NULL, only.preprocess=TRUE)
    }
    if (length(gopt))
        opt <- replace(opt, names(gopt), gopt)
    if (opt$RDataSuffix != suffix) {
        if (fix) {
            cat("Changing RDataSuffix saved in options to '", suffix, "'")
            if (dryRun)
                opt$RDataSuffix <- suffix
            else
                opt <- track.options(values=list(RDataSuffix=suffix), trackingEnv=trackingEnv)
            if (opt$RDataSuffix != suffix)
                stop("Strange: was unable to change stored option RDataSuffix to '", suffix, "'")
        } else {
            stop("Deduced RDataSuffix '", suffix, "' is different to suffix '",
                 opt$RDataSuffix,
                 "' stored in trackingEnv.  Do either track.rebuild(..., RDataSuffix='",
                 opt$RDataSuffix, "') to use the stored suffix, ",
                 "or track.rebuild(..., fix=TRUE) to change the stored suffix to '", suffix, "'")
        }
    }
    ##
    ## Now we have the options and directory and know if the db is active
    ##
    ## Read the list of actual objects in the env (if linked)
    if (activeTracking) {
        trackingEnvObjs <- ls(envir=trackingEnv, all.names=TRUE)
        if (length(trackingEnvObjs)) {
            i <- isReservedName(trackingEnvObjs)
            trackingEnvObjs <- trackingEnvObjs[!i]
        }
    } else {
        trackingEnvObjs <- character(0)
    }
    ## See if there are any files with a wrong suffix in the data directory
    ## If there are, then if fix==TRUE, try to rename these using the right
    ## suffix, or move them to the quarantine directory
    if (length(setdiff(gopt$RDataSuffixes, suffix))) {
        otherSuffixRegex <- paste("(", paste(setdiff(gopt$RDataSuffixes, suffix), collapse="|"), ")", sep="")
        files <- list.files(path=dataDir, pattern=paste(".*\\.", otherSuffixRegex, "$", sep=""), all.files=TRUE)
        ## If fix==TRUE, rename them or move them to quarantineDir
        if (length(files)) {
            cat("There are other ", length(files), " files in '", dataDir, "' that look like RDataFiles: ",
                paste("'", files[seq(len=min(2, length(files)))], "'", collapse=", ", sep=""),
                if (length(files)>2) "...", "\n", sep="")
            if (!fix || dryRun) {
                cat("Leaving these alone: call track.rebuild(..., fix=TRUE, dryRun=FALSE) to rename these with suffix '", suffix, "' or move them to quarantine\n", sep="")
            } else {
                cat("Attempting to rename these files to use suffix '", suffix, "'\n", sep="")
                for (f in files) {
                    f2 <- sub(paste("\\.", otherSuffixRegex, "$", sep=""), suffix, f)
                    if (file.exists(file.path(dataDir, f2))) {
                        if (file.exists(file.path(quarantineDir, f)))
                            cat("Cannot move file '", f, "' to quarantine: file of same name already exists there; continuing...\n", sep="")
                        else if (!file.rename(file.path(dataDir, f), file.path(quarantineDir, f)))
                            cat("Was unable to move file '", f, "' to quarantine; continuing...\n", sep="")
                        else
                            cat("Moved file '", f, "' to ", quarantineDir, "\n", sep="")
                    } else {
                        if (!file.rename(file.path(dataDir, f), file.path(dataDir, f2)))
                            cat("Was unable to renaming file '", f, "' to '", f2, "', continuing...\n", sep="")
                        else
                            cat("Renamed file '", f, "' to '", f2, "'\n", sep="")
                    }
                }
            }
        }
    }
    ## Use tmpenv for reading all saved objects
    tmpenv <- new.env(hash=TRUE, parent=emptyenv())
    ## Read the existing file map and summary, and get a list of the files
    ## and objects (in both the tracked env and tracking env), if the tracking
    ## db is active.
    ## Warn about objects that are in the file map or summary for which
    ## the object and the file don't exist (not much we can do about this).
    ## If fix==TRUE, delete these objects.

    ## Read the file map from filemap.txt
    ## open in binary mode so that we use just "\n" as a separator
    fileMapFound <- FALSE
    fileMap <- character(0)
    fileMapFullPath <- file.path(dataDir, "filemap.txt")
    if (file.exists(fileMapFullPath)) {
        open.res <- try(con <- file(fileMapFullPath, open="rb"), silent=TRUE)
        if (is(open.res, "try-error")) {
            cat("Cannot open '", abbrevWD(fileMapFullPath), "' for reading; ignoring and continuing... (error was: ",
                formatMsg(open.res), ")\n", sep="")
        } else {
            on.exit(close(con))
            fileData <- try(readLines(con=con, n=-1), silent=TRUE)
            if (is(fileData, "try-error")) {
                cat("Unable to read '", abbrevWD(fileMapFullPath), "'; ignoring and continuing... (error was: ",
                    formatMsg(fileData), ")\n", sep="")
            } else {
                ## Remove Windows line termination
                fileData <- gsub("\r", "", fileData)
                i <- regexpr(":", fileData, fixed=TRUE)
                if (any(i < 1)) {
                    cat("file map contains invalid data (need a ':' in each line): \"", abbrevWD(fileMapFullPath), "\": discarding this info and continuing...\n", sep="")
                } else {
                    fileMap <- substring(fileData, 1, i-1)
                    names(fileMap) <- substring(fileData, i+1)
                    fileMapFound <- TRUE
                }
            }
        }
    } else {
        cat("File map file ('", abbrevWD(fileMapFullPath), "') does not exist; continuing...\n", sep="")
    }

    if (activeTracking) {
        if (exists(".trackingFileMap", envir=trackingEnv, inherits=FALSE)) {
            ## See if the file map obj in the tracking env has other info...
            envFileMap <- get(".trackingFileMap", envir=trackingEnv, inherits=FALSE)
            if (!is.character(envFileMap)) {
                cat("Contents of .trackingFileMap in tracking env are not character data; discarding\n");
                envFileMap <- NULL
            } else if (is.null(names(envFileMap))) {
                cat("Contents of .trackingFileMap in tracking env do not have file names; discarding\n");
                envFileMap <- NULL
            } else if (any(names(envFileMap)=="")) {
                cat("Contents of .trackingFileMap in tracking env have some empty file names; discarding\n");
                envFileMap <- NULL
            } else if (any(is.na(names(envFileMap)))) {
                cat("Contents of .trackingFileMap in tracking env have some NA file names; discarding\n");
                envFileMap <- NULL
            } else if (any(envFileMap=="")) {
                cat("Contents of .trackingFileMap in tracking env have some empty obj names; discarding\n");
                envFileMap <- NULL
            } else if (any(is.na(envFileMap))) {
                cat("Contents of .trackingFileMap in tracking env have some NA obj names; discarding\n");
                envFileMap <- NULL
            }
            if (!is.null(envFileMap)) {
                if (length(fileMap)==0) {
                    cat("Found a fileMap in the environment, but not on disk\n")
                    fileMap <- envFileMap
                    fileMapFound <- TRUE
                } else {
                    cat("Found a fileMap in the environment and on disk\n")
                    i <- intersect(names(fileMap), names(envFileMap))
                    if (length(i) && any(fileMap[i] != envFileMap[i])) {
                        cat("Have conflicting mappings in environment and disk file map:\n")
                        print(cbind(environment=envFileMap[i], disk=fileMap[i]))
                        if (level != "all")
                            stop("can only proceed with level='all'")
                    }
                }
            } else {
                envFileMap <- character(0)
            }
        } else {
            if (fileMapFound)
                cat("Found a fileMap in the environment, but not on disk\n")
        }
    }
    ## Can now just use 'fileMap' and ignore 'envFileMap'

    ## Read the summary, into summary and envSummary
    ## Will keep both of these because and get most
    ## recent info for each object one at a time.
    summarySkel <- summaryRow(name="", opt=opt)
    envSummary <- summarySkel[0,]
    if (activeTracking) {
        if (exists(".trackingSummary", envir=trackingEnv, inherits=FALSE)) {
            envSummary <- getObjSummary(trackingEnv, opt=opt)
            if (!is.data.frame(envSummary)) {
                cat("Tracking summary from tracking environment is not a data frame - discarding it.\n")
                envSummary <- NULL
            }
            if (!identical(all.equal(names(envSummary), names(summarySkel)), TRUE)) {
                ## Get the summary columns in the right order
                if (identical(all.equal(sort(names(envSummary)), sort(names(summarySkel))), TRUE)) {
                    envSummary <- envSummary[,names(summarySkel),drop=FALSE]
                } else {
                    cat("Tracking summary from tracking environment has different column names - discarding it.\n")
                    envSummary <- summarySkel[0,]
                }
            }
            if (!is.null(envSummary))
                cat("Using tracking summary from tracking environment.\n")
        }
    }
    if (file.exists(file.path(dataDir, paste(".trackingSummary", suffix, sep=".")))) {
        load.res <- try(load(file=file.path(dataDir, paste(".trackingSummary", suffix, sep=".")), envir=tmpenv), silent=TRUE)
        if (is(load.res, "try-error")) {
            cat("Cannot load '", abbrevWD(file.path(dataDir, paste(".trackingSummary", opt$RDataSuffix, sep="."))),
                "' -- rebuilding... (error was: ",
                formatMsg(load.res), ")\n")
        } else if(!is.element(".trackingSummary", load.res)) {
            cat("Strange: '", abbrevWD(file.path(dataDir, paste(".trackingSummary", opt$RDataSuffix, sep="."))),
                " does not contain a '.trackingSummary' object -- rebuilding it ...\n")
        } else {
            fileSummary <- getObjSummary(tmpenv, opt=opt)
            if (!is.data.frame(fileSummary)) {
                cat("Tracking summary read from file is not a data frame - ignoring it.\n")
                fileSummary <- NULL
            }
            if (!identical(all.equal(names(fileSummary), names(summarySkel)), TRUE)) {
                if (identical(all.equal(sort(names(fileSummary)), sort(names(summarySkel))), TRUE)) {
                    fileSummary <- fileSummary[,names(summarySkel),drop=FALSE]
                } else {
                    cat("Tracking summary read from file has different column names - ignoring it.\n")
                    fileSummary <- NULL
                }
            }
        }
        if (is.null(fileSummary))
            fileSummary <- summarySkel[0,]
        else
            cat("Using tracking summary read from file.\n")

        objs <- ls(all.names=TRUE, envir=tmpenv)
        if (length(objs))
            remove(list=objs, envir=tmpenv)
    } else {
        fileSummary <- summarySkel[0,]
    }
    ## After now, use both 'fileSummary' and 'envSummary'

    ## Read the list of actual files in the db.
    dbFiles <- list.files(dataDir, paste(".*\\.", suffix, sep=""), all.files=TRUE)
    names(dbFiles) <- sub(paste("\\.", suffix, "$", sep=""), "", dbFiles)
    dbFiles <- dbFiles[!is.element(names(dbFiles), c(".trackingSummary", ".trackingOptions"))]
    ## The contents of dbFiles are the file names with suffix, the names
    ## are the file names without suffixes.

    ## Work out which files we need to re-read (if level='missing'):
    ##   (1) files that are not mentioned in fileMap
    ##   (2) files for objects that are not in the summary
    ## If level='missing' and there are files in the fileMap
    ## that are not in the directory, and there is no corresponding
    ## object in the trackingEnv, stop with an error, because
    ## continuing will result in the meta-information for these
    ## objects being lost.
    ## Re-read all files if level='all'
    if (length(setdiff(fileMap, names(dbFiles)))) {
        i <- match(setdiff(fileMap, names(dbFiles)), fileMap)
        objs <- names(fileMap)[i]
        ## if these objects are in the tracking env, proceed without comment
        objs <- setdiff(objs, trackingEnvObjs)
        if (length(objs)) {
            i <- match(objs, names(fileMap))
            cat("Cannot find some objects in the fileMap: no file on disk, and object is not in the tracking environment:\n")
            print(data.frame(filename=paste(fileMap[i], suffix, sep="."), row.names=names(fileMap)[i]))
            if (level=="missing")
                stop("can only proceed with level='all' -- and in that case the meta info for these lost objects will be lost")
        }
    }

    ## See what information from the file map and summary we can completely
    ## reuse (and not have to look at the object)
    summaryObjNames <- unique(c(row.names(envSummary), row.names(fileSummary)))
    # reuseFileMap <- NULL # don't need this object
    reuseSummary <- NULL
    if (level=="missing") {
        ## Just read files that are not in the filemap
        filesToRead <- setdiff(names(dbFiles), fileMap)
        if (length(filesToRead)) {
            cat("Some data files are not in the file map: ", paste("'", filesToRead, ".", suffix, "'", sep="", collapse=", "), sep="", "\n")
        }
        if (any(i <- !is.element(names(fileMap), summaryObjNames))) {
            cat("Some objects in the file map are not in the summary:\n")
            cat(paste(format(c("Object", names(fileMap)[i])), " ", format(c("File", paste(fileMap[i], suffix, sep="."))), "\n", sep=""), sep="")
            filesToRead <- c(filesToRead, fileMap[i])
        }
    } else {
        filesToRead <- unique(c(names(dbFiles), fileMap))
        reuseObjs <- character(0)
    }
    ## Reuse data for objects in the fileMap that are in the summary
    ## and have objects or files.
    ## First look at objects in the env
    objNames <- intersect(names(fileMap), row.names(envSummary))
    reuseSummary <- envSummary[objNames,,drop=FALSE]
    objNames <- setdiff(intersect(names(fileMap), row.names(fileSummary)), objNames)
    if (length(objNames))
        reuseSummary <- rbind(reuseSummary, fileSummary[objNames,,drop=FALSE])
    reuseObjs <- row.names(reuseSummary)
    # reuseFileMap <- fileMap[reuseObjs]

    ##
    ## Read each of the files in filesToRead (or its object in the env)
    ## Note that an object mentioned here may not exist on disk -- it
    ## may only exist in the environment (if was cached and not yet saved
    ## when track.rebuild() was called.)  In general, we assume that a
    ## cached version of an object is the one to use.
    ##
    newFileMap <- character(length(filesToRead))
    objsFromEnv <- character(0)
    sumRowsReused <- 0
    newSummaryRows <- vector("list", length(filesToRead))
    for (objFileBase in filesToRead) {
        objFile <- paste(objFileBase, ".", suffix, sep="")
        objName <- NA
        objEnv <- NULL
        if (!is.na(i <- match(objFileBase, fileMap)))
            objName <- names(fileMap)[i]
        if (activeTracking && is.element(objName, trackingEnvObjs)) {
            if (verbose > 1)
                cat("Using cached version of object '", objName, "'")
            ## Find the object if it exists, and use it in preference
            ## to the file copy (which could be different and possibly
            ## out of date.)
            objEnv <- trackingEnv
            objNames <- objName
            objsFromEnv <- c(objsFromEnv, objName)
        } else if (!file.exists(file.path(dataDir, objFile))) {
            cat("Object '", (if (!is.na(objName)) objName  else "(unknown)") , "'",
                " is not stored on file ('", abbrevWD(file.path(dataDir, objFile)), "')",
                " and is not cached -- forgetting about it\n", sep="")
            i <- match(objName, names(fileMap))
            if (!is.na(i))
                fileMap <- fileMap[-i]
            ## Remove the active binding
            if (!dryRun && activeTracking && exists(objName, envir=envir, inherits=FALSE))
                rm(list=objName, envir=envir)
            if (is.element(objName, rownames(reuseSummary)))
                reuseSummary <- reuseSummary[-match(objName, rownames(reuseSummary)), , drop=FALSE]
        } else {
            if (verbose > 1)
                cat("Loading file '", abbrevWD(file.path(dataDir, objFile)), "' for ",
                    (if (!is.na(objName)) paste("object '", objName, "'", sep="") else "unknown object"),
                    " ... ", sep="")
            load.res <- try(load(file=file.path(dataDir, objFile), envir=tmpenv), silent=TRUE)
            if (is(load.res, "try-error")) {
                cat(": loading error\n")
            } else {
                cat("loaded ", paste("'", load.res, "'", sep="", collapse=", "), "\n", sep="")
            }
            ok <- TRUE
            fixable <- TRUE
            if (dryRun) {
                moving.msg <- "would move"
                rewriting.msg <- "rewrite"
            } else {
                moving.msg <- "moving"
                rewriting.msg <- "rewriting"
            }
            if (is(load.res, "try-error")) {
                cat(objFile, " cannot be loaded (error was '", formatMsg(load.res), "')",
                    if (fix) paste(",", moving.msg, "to", abbrevWD(quarantineDir)), "\n", sep="")
                ok <- FALSE
                fixable <- FALSE
            } else if (length(load.res)!=1) {
                cat(objFile, " cannot be used as is because it contains multiple objects",
                    if (fix) paste(",", moving.msg, "to", abbrevWD(quarantineDir), "and", rewriting.msg, "objects\n"), sep="")
                ok <- FALSE
            } else if (is.element(load.res, names(newFileMap))) {
                cat(objFile, " cannot be used because it contains an object with the same name",
                    " as an object in another newly read file ('",
                    load.res, "', also in ", newFileMap[load.res], ")",
                    if (fix) paste(",", moving.msg, "to", abbrevWD(quarantineDir)), "\n", sep="")
                ok <- FALSE
            } else if (is.element(load.res, names(fileMap))) {
                # here we have an object that was mentioned in the filemap
                if (is.na(objName)) {
                    # this file wasn't in the filemap
                    cat(objFile, " cannot be used because it contains an object with the same name",
                        " as an object in a different file mentioned in the existing filemap ('",
                        load.res, "', also in ", fileMap[load.res], ")",
                        if (fix) paste(",", moving.msg, "to", abbrevWD(quarantineDir)), "\n", sep="")
                    ok <- FALSE
                } else if (objName != load.res) {
                    # The filename was in the filemap, but the object is different.
                    # This is a bizarre situation -- it looks like the save file has
                    # been switched to contain a different object.
                    # May could recover here by re", moving.msg, " reference to the old object...,
                    # but for the moment, code this as a move-to-quarantine situation.
                    cat(objFile, " cannot be used because it was expected to contain the object named '",
                        load.res, "'",
                        if (fix) paste(",", moving.msg, "to", abbrevWD(quarantineDir)), "\n", sep="")
                    ok <- FALSE
                }
            } else if (isSimpleName(objFileBase) && objFileBase!=load.res) {
                cat(objFile, " cannot be used as is because the file has a 'simple' name different to the object name",
                    if (fix) paste(",", moving.msg, "to", abbrevWD(quarantineDir), "and", rewriting.msg, "the object"), "\n", sep="")
                ok <- FALSE
            } else if (!isSimpleName(load.res) && substring(objFileBase, 1, 1) != "_") {
                cat(objFile, " cannot be used as is because objects that do not have a 'simple' name must be stored in a file beginning with '_'",
                    if (fix) paste(",", moving.msg, "to", abbrevWD(quarantineDir), "and", rewriting.msg, "the object"), "\n", sep="")
                ok <- FALSE
            } else if (isSimpleName(load.res) && load.res!=objFileBase && substring(objFileBase, 1, 1) != "_") {
                cat(objFile, " cannot be used as is because objects that do not have a 'simple' name must be stored in a file beginning with '_'",
                    if (fix) paste(",", moving.msg, "to", abbrevWD(quarantineDir), "and", rewriting.msg, "the object"), "\n", sep="")
                ok <- FALSE
            }
            if (ok) {
                objEnv <- tmpenv
                objNames <- load.res
                ## Add to the new filemap if it wasn't in the old filemap
                if (is.na(objName))
                    newFileMap <- setNamedElt(newFileMap, load.res, objFileBase)
                ## Call the load callback for this object
            } else if (fix) {
                if (!dryRun) {
                    if (!dir.exists(quarantineDir))
                        dir.create(quarantineDir)
                    if (!file.rename(file.path(dataDir, objFile), file.path(quarantineDir, objFile)))
                        cat("Unable to move '", objFile, "' to quarantine\n", sep="")
                }
                if (fixable) {
                    ## Write the object out into file(s) with legal names.
                    ## Include names(dbFiles) to make sure that the new names
                    ## don't conflict with files in the dir
                    newFileBase <- makeObjFileName(load.res, list(fileMap, newFileMap, names(dbFiles)))
                    savedOk <- logical(length(load.res))
                    for (i in seq(along=load.res)) {
                        if (is.element(load.res[i], names(newFileMap))) {
                            cat("Ignoring object with already-used name ('", load.res, "' in '", objFile, "'\n", sep="")
                        } else {
                            newFile <- file.path(dataDir, paste(newFileBase[i], suffix, sep="."))
                            cat(if (dryRun) "Would save" else "Saving",
                                " object '", load.res[i], "' to '", abbrevWD(newFile), "'\n", sep="")
                            if (!dryRun) {
                                save.res <- try(save(list=load.res[i], file=newFile, envir=tmpenv,
                                                     compress=opt$compress, compression_level=opt$compression_level), silent=TRUE)
                                if (is(save.res, "try-error")) {
                                    cat("Could not save obj '", load.res[i], "' in file '", abbrevWD(newFile), "' (error was '",
                                        formatMsg(save.res), "')\n", sep="")
                                } else {
                                    savedOk[i] <- TRUE
                                }
                            } else {
                                # dryRun: pretend that it was saved OK
                                savedOk[i] <- TRUE
                            }
                        }
                    }
                    newFileMap <- setNamedElt(newFileMap, load.res[savedOk], newFileBase[savedOk])
                    objEnv <- tmpenv
                    objNames <- load.res[savedOk]
                }
            } else {
                stop("Cannot proceed further -- rerun with fix=TRUE")
            }
        }
        if (!is.null(objEnv)) {
            ## objEnv could be the tracking env or tmpenv
            ## if objEnv is NULL, means we could not find the object -- skip it
            for (o in objNames) {
                obj <- get(o, envir=objEnv, inherits=FALSE)
                sumRow <- chooseBestSummaryRow(o)
                if (is.null(sumRow)) {
                    ## try to rebuild this summary row
                    if (use.file.times)
                        sumRow <- summaryRow(o, opt=opt, obj=obj, file=file.path(dataDir, objFile), accessed=FALSE)
                    else
                        sumRow <- summaryRow(o, opt=opt, obj=obj, accessed=FALSE)
                } else {
                    ## update the summary row (but don't change times)
                    sumRow <- summaryRow(o, opt=opt, sumRow=sumRow, obj=obj, accessed=FALSE)
                    sumRowsReused <- sumRowsReused + 1
                }
                newSummaryRows <- setNamedElt(newSummaryRows, o, sumRow)
            }
        }
        ## Clean up tmpenv for the next iteration
        objs <- ls(all.names=TRUE, envir=tmpenv)
        if (length(objs))
            remove(list=objs, envir=tmpenv)
    }

    ## If there are any objects in the tracking env that we haven't looked at,
    ## look at those now.  These will need entries in the fileMap and probably
    ## the summary too (otherwise we would have looked at them already.)
    ## Add entries to .trackingUnsaved for these vars (don't save them to files)
    unsaved <- NULL
    if (activeTracking) {
        objsToDo <- setdiff(trackingEnvObjs, c(objsFromEnv, reuseObjs))
        objsToDo <- objsToDo[!isReservedName(objsToDo)]
        if (length(objsToDo)) {
            for (objName in objsToDo) {
                if (verbose > 1)
                    cat("Looking at object '", objName, "' in tracking environment\n", sep="")
                if (is.element(objName, names(fileMap)))
                    stop("should not happen: didn't expect to find '", objName, "' in fileMap")
                if (is.element(objName, names(newFileMap)))
                    stop("should not happen: didn't expect to find '", objName, "' in newFileMap")
                newFileBase <- makeObjFileName(objName, list(fileMap, newFileMap))
                newFileMap <- setNamedElt(newFileMap, objName, newFileBase)
                obj <- get(objName, inherits=FALSE, envir=trackingEnv)
                sumRow <- chooseBestSummaryRow(objName)
                ## update or recreate the summary row (sumRow could be NULL)
                ## we have no file info to get times from
                sumRow <- summaryRow(objName, opt=opt, sumRow=sumRow, obj=obj)
                newSummaryRows[[objName]] <- sumRow
            }
            unsaved <- mget(".trackingUnsaved", envir=trackingEnv, inherits=FALSE, ifnotfound=list(character(0)))[[1]]
            unsaved <- sort(unique(c(unsaved, objsToDo)))
            assign(".trackingUnsaved", unsaved, envir=trackingEnv)
        }
    }
    newSummary <- do.call("rbind", newSummaryRows)
    ## Get just the elements of newFileMap that we actually used.
    if (is.null(names(newFileMap)))
        newFileMap <- character(0)
    else
        newFileMap <- newFileMap[!is.na(names(newFileMap))]

    ## Combine the old and new info.
    if (verbose>1)
        cat("New file map has", length(newFileMap), "newly constructed entries and uses", length(fileMap), "old entries\n")
    newFileMap <- c(newFileMap, fileMap)
    if (length(newFileMap))
        newFileMap <- newFileMap[order(names(newFileMap))]
    if ((!is.null(newSummary) && nrow(newSummary)) || (!is.null(reuseSummary) && nrow(reuseSummary))) {
        if (verbose>1)
            cat("New summary has",
                (if (is.null(newSummary)) 0 else nrow(newSummary)) - sumRowsReused,
                "newly constructed entries and uses",  sumRowsReused, "recovered and",
                (if (is.null(reuseSummary)) 0 else nrow(reuseSummary)),
                "old entries\n")
        newSummary <- rbind(reuseSummary, newSummary[!(rownames(newSummary) %in% rownames(reuseSummary)), , drop=FALSE])
    } else {
        if (verbose>1)
            cat("New summary has zero entries\n")
        newSummary <- summarySkel[0,]
    }
    if (nrow(newSummary))
        newSummary <- newSummary[order(row.names(newSummary)),,drop=FALSE]

    ## See if any objects are missing their active bindings
    masked <- character(0)
    if (activeTracking) {
        ## Get all the objects in the environment
        visibleObjs <- ls(envir=envir, all.names=TRUE)
        ## We're only interested in ones with the same names as tracked vars
        visibleObjs <- intersect(visibleObjs, rownames(newSummary))
        missingBindings <- setdiff(rownames(newSummary), visibleObjs)
        if (length(missingBindings)) {
            cat(if (dryRun) "Would rebuild " else "Rebuilding ",
                length(missingBindings), " missing active bindings: ",
                paste(missingBindings, collapse=", "), "\n", sep="")
            if (!dryRun) {
                for (objName in missingBindings) {
                    f <- createBindingClosure(objName, trackingEnv)
                    makeActiveBinding(objName, env=envir, fun=f)
                }
            }
        }
        visibleObjsAB <- sapply(visibleObjs, bindingIsActive, envir)
        if (any(!visibleObjsAB)) {
            cat(sum(!visibleObjsAB), "tracked vars are masked by non-active bindings -- to fix, please manually remove these vars from the tracked environment before re-running track.rebuild()\n")
            cat("The masked vars are: ", paste(visibleObjs[!visibleObjsAB], collapse=", "), "\n", sep="")
            masked <- visibleObjs[!visibleObjsAB]
        }
    }

    res <- list(fileMap=newFileMap, summary=newSummary)
    if (length(unsaved))
        res$unsaved <- unsaved
    if (length(masked))
        res$masked <- masked
    if (dryRun) {
        cat("Run with dryRun=FALSE to actually make changes\n")
        return(invisible(res))
    } else {
        if (activeTracking) {
            if (verbose>1)
                cat("Assigning '.trackingFileMap' and '.trackingSummary' in tracking environment.\n")
            if (is(assign.res <- try(assign(".trackingFileMap", newFileMap, envir=trackingEnv), silent=TRUE), "try-error"))
                warning("failed to assign '.trackingFileMap' in ", envname(trackingEnv),
                        " (error was '", formatMsg(res), "')")
            assign.res <- try(assign(".trackingSummary", newSummary, envir=trackingEnv), silent=TRUE)
            if (is(assign.res, "try-error"))
                warning("unable to assign .trackingSummary back to tracking env on ", envname(envir),
                        " (error was '", formatMsg(assign.res), "')")
            else
                assign(".trackingSummaryChanged", TRUE, envir=trackingEnv)
        }
        if (verbose>1)
            cat("Writing file map to file.\n")
        writeFileMapFile(newFileMap, trackingEnv, dataDir, FALSE)
        file <- file.path(getDataDir(dir), paste(".trackingSummary", opt$RDataSuffix, sep="."))
        if (verbose>1)
            cat("Saving object summary to file.\n")
        if (activeTracking) {
            save.res <- saveObjSummary(trackingEnv, opt, getDataDir(dir))
        } else {
            ## Need to have newSummary in a variable named '.trackingSummary' to use save()
            assign(".trackingSummary", newSummary, envir=tmpenv)
            save.res <- saveObjSummary(trackingEnv, envir=tmpenv, opt=opt, dataDir=getDataDir(dir))
            # save.res <- try(save(list=".trackingSummary", file=file, envir=tmpenv, compress=FALSE), silent=TRUE)
        }
        if (is(save.res, "try-error"))
            warning("unable to save .trackingSummary to ", dir, " (error was '", formatMsg(save.res), "')")
        else if (activeTracking && !is(assign.res, "try-error"))
            assign(".trackingSummaryChanged", FALSE, envir=trackingEnv)
        return(invisible(res))
    }
}

