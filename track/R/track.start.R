track.start <- function(dir="rdatadir", pos=1, envir=as.environment(pos),
                        create=TRUE, clobber=c("no", "files", "variables", "vars", "var"),
                        discardMissing=FALSE,
                        cache=NULL, cachePolicy=NULL,
                        options=NULL, RDataSuffix=NULL, auto=NULL,
                        readonly=FALSE, lockEnv=FALSE, check.Last=TRUE,
                        autoCheckSize=1e6,
                        verbose=TRUE) {
    ## Start tracking the specified environment to a directory
    clobber <- match.arg(clobber)
    if (clobber=="vars" || clobber=="var") clobber <- "variables"
    if (env.is.tracked(envir))
        stop("env ", envname(envir), " is already tracked by dir '",
             get(".trackingDir", envir=getTrackingEnv(envir), inherits=FALSE), "'")
    dir.orig <- dir
    dir <- getAbsolutePath(dir)
    ## Working out the options values to use is a little tricky.
    ## This is the priority:
    ## (1) values in the 'options' argument ('cache' override options$cache)
    ## (2) values in a saved .trackOptions object (if it exists)
    ## (3) getOption("global.track.options") (user modifiable)
    ## (4) standard package defaults (not user modifiable)
    ##
    ## Create the tracking env, but don't assign it to envir until all is OK
    trackingEnv <- new.env(hash=TRUE, parent=emptyenv())
    assign(".trackingDir", dir, envir=trackingEnv)
    track.stop.finalizer <- function(trackingEnv) {
        ## Finalizer is difficult because it can be called long
        ## after a tracking environment has been disconnected.
        ## Had the following disabled checks in here because I thought I was
        ## seeing some problems with invalid calls to the finalizer, but they
        ## all turned out to be a result of the finalizer on the object being
        ## called long after it had stopped being used.
        if (exists(".trackingFinished", envir=trackingEnv, inherits=FALSE))
            return(NULL)
        if (!exists(".trackingEnv", envir=envir, inherits=FALSE)) {
            ## This used to happen under some circumstances when the finalizer is
            ## called after tracking has stopped, but the check for ".trackingFinished"
            ## fixed that.
            if (FALSE)
                cat("Bogus call to track.stop reg.finalizer for", envname(trackingEnv),
                    ": no .trackingEnv in", envname(envir), "\n")
            return(NULL)
        }
        if (!identical(trackingEnv, get(".trackingEnv", envir=envir, inherits=FALSE))) {
            ## This can happen in cases where a tracking env is partially
            ## set up and then discarded, as when there are variable name
            ## conflicts
            if (FALSE)
                cat("Bogus call to track.stop reg.finalizer for", envname(trackingEnv),
                    ": .trackingEnv in", envname(envir), "is different:",
                    envname(get(".trackingEnv", envir=envir, inherits=FALSE)), "\n")
            return(NULL)
        }
        ## cat("Valid call to track.stop reg.finalizer for", envname(trackingEnv),
        ##         "on", envname(envir), "\n")
        ## if (interactive()) browser()
        ## Seems to sometimes be called when not appropriate, so don't
        ## do anything drastic -- just flush.
        ## Example is the track.load() commands in the examples in track.status.Rd
        track.flush(envir=envir)
    }
    if (!readonly) {
        ## cat("Installing reg.finalizer for trackingEnv", envname(trackingEnv), "on", envname(envir), "\n")
        reg.finalizer(trackingEnv, track.stop.finalizer, onexit=TRUE)
    }
    ## Set up a finalizer to be run when the tracking env is garbage collected
    ## OR at the end of the session
    dataDir <- getDataDir(dir)
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
    if (is.null(auto))
        if (length(gopt$autoTrack))
            auto <- gopt$autoTrack
        else
            auto <- TRUE
    optionsPath <- NULL
    if (dir.exists(file.path(dataDir))) {
        ## Try to work out the suffix being used
        ## First look for .trackingOptions file
        suffix <- NULL
        x <- list.files(path=dataDir, pattern=paste("^\\.trackingOptions\\.", suffixRegExp, "$", sep=""), all.files=TRUE)
        if (length(x)>1)
            stop("have multiple options files in '", dataDir, "': ", paste(x, collapse=", "))
        if (length(x)==1) {
            optionsPath <- file.path(dataDir, x)
            suffix <- sub(".*\\.", "", x)
        }
        if (is.null(suffix)) {
            ## next look for .trackingSummary file
            x <- list.files(path=dataDir, pattern=paste("^\\.trackingSummary\\.", suffixRegExp, "$", sep=""), all.files=TRUE)
            if (length(x)>1)
                stop("have multiple summary files in '", dataDir, "': ", paste(x, collapse=", "))
            if (length(x)==1)
                suffix <- sub(".*\\.", "", x)
        }
        if (is.null(suffix)) {
            ## next look for any files with possible RData suffix
            x <- list.files(path=dataDir, pattern=paste("^.*\\.", suffixRegExp, "$", sep=""), all.files=TRUE)
            if (length(x)>0) {
                suffix <- unique(sub(".*\\.", "", x))
                if (length(suffix)>1)
                    stop("have files with multiple RData suffixes in '", dataDir, "': ", paste(x, collapse=", "))
            } else {
                if (!is.null(RDataSuffix))
                    suffix <- RDataSuffix
                else
                    suffix <- gopt$RDataSuffixes[1]
            }
        }
        if (verbose)
            cat("Tracking ", envname(envir), if (readonly) " (readonly)" else " (writable)",
                " using existing directory '", dir.orig, "'\n", sep="")
    } else {
        if (is.null(RDataSuffix))
            suffix <- gopt$RDataSuffixes[1]
        else
            suffix <- RDataSuffix
        if (verbose)
            cat("Tracking ", envname(envir), if (readonly) " (readonly)" else " (writable)",
                " using new directory '", dir.orig, "'\n", sep="")
    }
    if (!is.element(suffix, gopt$RDataSuffixes))
        stop("internal error: ended up with an illegal suffix?? (", suffix, ")")
    if (!is.null(RDataSuffix) && RDataSuffix != suffix)
        stop("suffix in use '", suffix, "' differs from supplied RDataSuffix ('", RDataSuffix, "')")

    ## Preprocess the options (standardize, get defaults).
    if (!is.null(cache)) {
        if (is.null(options))
            options <- list()
        if (is.logical(cache) && length(cache)==1 && !is.na(cache))
            options$cache <- cache
        else
            stop("'cache' argument must be TRUE or FALSE")
    }
    if (!is.null(cachePolicy)) {
        if (is.null(options))
            options <- list()
        if (is.element(cachePolicy, c("eotPurge", "none")))
            options$cachePolicy <- cachePolicy
        else
            stop("'cachePolicy' argument must be 'eotPurge' or 'none'")
    }
    old.options <- list()
    if (!is.null(optionsPath)) {
        tmpenv <- new.env(parent=emptyenv())
        load.res <- try(load(optionsPath, envir=tmpenv), silent=TRUE)
        if (is(load.res, "try-error"))
            warning(optionsPath, " cannot be loaded -- ignoring it; for recovery see ?track.rebuild (",
                 as.character(load.res), ")")
        if (length(load.res)!=1 || load.res != ".trackingOptions") {
            warning(optionsPath, " does not contain just '.trackingOptions' -- ignoring it; for recovery see ?track.rebuild")
        } else {
            ## .trackingOptions has to exist because we just loaded it
            old.options <- get(".trackingOptions", envir=tmpenv, inherits=FALSE)
            if (!is.list(old.options)) {
                warning("'.trackingOptions' from ", optionsPath, " is not a list -- ignoring it; see ?track.rebuild")
                old.options <- list()
            }
        }
    }
    if (length(old.options)==0) {
        ## Couldn't read any options from the file, initialize them
        ## from remaining gopt.
        gopt$autoTrack <- NULL
        gopt$RDataSuffixes <- NULL
        old.options <- gopt
    }
    opt <- track.options(values=options, envir=NULL, only.preprocess=TRUE, old.options=old.options)
    if (!is.null(readonly))
        opt$readonly <- readonly
    assign(".trackingOptions", opt, envir=trackingEnv)
    fileMapPath <- file.path(dataDir, "filemap.txt")
    fileMapChanged <- FALSE
    ## Create a default empty objSummary -- an existing one will replace this
    objSummary <- summaryRow(name="", opt=opt)[0,]
    if (!dir.exists(file.path(dataDir))) {
        if (!create)
            stop("dir \"", dataDir, "\" does not exist (supply create=TRUE to create it)")
        res <- dir.create(file.path(dataDir), recursive=TRUE)
        if (is(res, "try-error"))
            stop("could not creating tracking dir '", dataDir, "': ", res)
        if (!dir.exists(dataDir))
            stop("failed to create tracking dir '", dataDir, "'")
        fileMap <- character(0)
        assign(".trackingFileMap", fileMap, envir=trackingEnv)
        fileMapChanged <- TRUE
        assign(".trackingSummary", objSummary, envir=trackingEnv)
        if (opt$debug >= 2)
            cat('track.start: created new tracking dir and assigned fileMap and objSummary in trackingEnv\n')
    } else {
        if (opt$debug >= 2)
            cat('track.start: attempting to read summary and fileMap from existing tracking db\n')
        ## Try to read the summary first, because if there is a problem with fileMap,
        ## we may want to erase the summary (though the code doesn't currently
        ## do that.)
        objSummaryFromFile <- loadObjSummary(trackingEnv, opt, stop.on.not.exists=FALSE)
        if (!is.null(objSummaryFromFile)) {
            objSummary <- objSummaryFromFile
            assign(".trackingSummary", objSummary, envir=trackingEnv)
        }
        ## Need to confirm that there are no clashes between variables already
        ## stored in the tracking dir and variables in envir
        if (file.exists(fileMapPath)) {
            fileMap <- readFileMapFile(trackingEnv, dataDir, TRUE)
            if (length(fileMap)) {
                fileExists <- file.exists(file.path(dataDir, paste(fileMap, sep=".", opt$RDataSuffix)))
                if (any(!fileExists)) {
                    if (discardMissing) {
                        if (verbose)
                            cat('Discarding info about objects with missing save files: ',
                                 paste(names(fileMap)[!fileExists], collapse=", "), "\n", sep="")
                        if (any(names(fileMap)[!fileExists] %in% rownames(objSummary)))
                            objSummary <- objSummary[(rownames(objSummary) %in% names(fileMap)[fileExists]), , drop=FALSE]
                        fileMap <- fileMap[fileExists]
                        fileMapChanged <- TRUE
                    } else {
                        warning("missing files for some variables in the fileMap (supply discardMissing=TRUE or remove or assign variables to repair): ",
                                paste(names(fileMap)[!fileExists], collapse=", "))
                    }
                }
            }
            alreadyExists <- logical(0)
            if (length(fileMap))
                alreadyExists <- sapply(names(fileMap), exists, envir=envir, inherits=FALSE)
            alreadyExists <- names(fileMap)[alreadyExists]
            ## opt$clobberVars contains a vector of var names that it is
            ## always OK to clobber.
            ## opt$clobberVars is processed before the clobber= argument
            clobberFirst <- intersect(opt$clobberVars, alreadyExists)
            if (length(clobberFirst)) {
                alreadyExists <- setdiff(alreadyExists, clobberFirst)
                remove(list=clobberFirst, envir=envir)
            }
            i <- FALSE
            for (re in opt$autoTrackExcludePattern)
                i <- i | grep(re, names(fileMap))
            if (any(i))
                warning("tracking db contains some vars that match the autoExclude pattern (this will be tracked, and continue to be tracked until removed): ", names(fileMap)[i])
            i <- isReservedName(names(fileMap))
            if (any(i))
                warning("tracking db contains some vars that have reserved names (this shouldn't happen, and may affect the correct operation of tracking): ", names(fileMap)[i])
            if (length(alreadyExists)) {
                if (clobber=="no") {
                    # If the objects and files are small, see if they are the same
                    tmpenv <- new.env(parent=emptyenv())
                    objSizes <- sapply(alreadyExists, function(v) object.size(get(v, envir=envir, inherits=FALSE)))
                    fileSizes <- file.info(file.path(dataDir, paste(fileMap, sep='.', opt$RDataSuffix)))$size
                    knowSame <- rep(FALSE, length(alreadyExists))
                    if (sum(objSizes) < autoCheckSize && sum(fileSizes, na.rm=TRUE) < autoCheckSize) {
                        knowSame <- sapply(alreadyExists, function(objName) {
                            objFile <- file.path(dataDir, paste(fileMap[objName], sep='.', opt$RDataSuffix))
                            load.res <- try(load(objFile, envir=tmpenv), silent=TRUE)
                            if (is(load.res, "try-error"))
                                stop("Failed to load R object ", objName, " from file ", objFile,
                                     " when checking whether existing R objects are same as those in the",
                                     " tracking db on the file system -- repair or delete file and try again;",
                                     " problem was: ", as.character(load.res))
                            if (length(load.res)!=1 || load.res != objName)
                                stop(objFile, " does not contain just '", objName, "' -- for recovery see ?track.rebuild")
                            objValueFile <- get(objName, envir=tmpenv, inherits=FALSE)
                            objValueEnv <- get(objName, envir=envir, inherits=FALSE)
                            rm(list=objName, envir=tmpenv, inherits=FALSE)
                            return(identical(objValueEnv, objValueFile))
                        })
                        if (all(knowSame)) {
                            # Must remove objs from the env to make way for the active bindings
                            # Todo: if some vars have special treatment such as alwaysCache,
                            # then transfer them into the tracking environment.
                            if (verbose)
                                warning('Have identical vars with same names in tracking db in "', dir, '" and in ',
                                        envname(envir), ': ', paste("'", alreadyExists, "'", sep='', collapse=', '))
                            remove(list=alreadyExists, envir=envir)
                        } else {
                            # If not all vars are the same, put ones that differ at the start of
                            # 'alreadyExists' to make for a more informative error message.
                            alreadyExists <- alreadyExists[order(knowSame, alreadyExists)]
                        }
                    }
                    if (!all(knowSame)) {
                        assign(".trackAlreadyExists", alreadyExists, envir=envir)
                        stop("cannot start tracking to dir \"", dir, "\" because it contains ",
                             length(alreadyExists), " vars that currently exist in ", envname(envir),
                             ", e.g.: ", paste("'", alreadyExists[seq(len=min(3,length(alreadyExists)))], "'", sep="", collapse=", "),
                             if (length(alreadyExists)>3) ", ...",
                             " (try track.start(..., clobber='files') or track.start(..., clobber='vars') to clobber one or the other")
                    }
                } else if (clobber=="files") {
                    if (opt$readonly) {
                        warning("will not clobber files corresponding to existing variables because readonly=TRUE: ", paste(alreadyExists, collapse=", "))
                    } else {
                        for (varName in alreadyExists) {
                            file <- fileMap[match(varName, names(fileMap))]
                            file.remove(file.path(dataDir, paste(file, opt$RDataSuffix, sep=".")))
                            value <- get(varName, envir=envir, inherits=FALSE)
                            setTrackedVar(varName, value, trackingEnv, opt=replace(opt, "maintainSummary", FALSE), file=file)
                            remove(list=varName, envir=envir, inherits=FALSE)
                        }
                    }
                } else if (clobber=="variables") {
                    remove(list=alreadyExists, envir=envir)
                }
            }
        } else {
            ## if there are in .rda files in this directory, need to rebuild the fileMap
            files <- list.files(path=file.path(dataDir), pattern=paste(".*\\.", opt$RDataSuffix, "$", sep=""), all.files=TRUE)
            files <- setdiff(files, paste(".trackingSummary", opt$RDataSuffix, sep="."))
            if (length(files)==0) {
                ## No files: start with a new fileMap and objSummary
                fileMap <- character(0)
                fileMapChanged <- TRUE
                assign(".trackingFileMap", fileMap, envir=trackingEnv)
            } else {
                stop("tracking dir \"", dir, "\" has some data files in it, but has no filemap.txt file -- use track.rebuild() to fix it")
            }
        }
        if (nrow(objSummary)) {
            ## update the prior-reads and prior-writes, and existing session counts
            objSummary[,"ES"] <- objSummary[,"ES"] + 1
            objSummary[,"PA"] <- objSummary[,"PA"] + objSummary[,"SA"]
            objSummary[,"PW"] <- objSummary[,"PW"] + objSummary[,"SW"]
            objSummary[,"SA"] <- 0
            objSummary[,"SW"] <- 0
        }
    }
    ## Do we need to save the fileMap ? (only if changed)
    if (fileMapChanged && !opt$readonly)
        writeFileMapFile(fileMap, trackingEnv, dataDir, TRUE)
    ## We always need to save the summary...
    assign(".trackingSummary", objSummary, envir=trackingEnv)
    assign(".trackingSummaryChanged", TRUE, envir=trackingEnv)
    if (dir!=dataDir) {
        ## Only use a "DESCRIPTION" file when we use a 'data' subdirectory
        ## in the tracking dir.
        dfn <- file.path(dir, "DESCRIPTION")
        if (!file.exists(dfn)) {
            write.res <- try(cat(track.package.desc(basename(dir)), "\n", sep="\n", file=dfn), silent=TRUE)
            if (is(write.res, "try-error"))
                warning("had problem writing ", dfn, " (", as.character(write.res), ")")
        }
    }
    ## create bindings for the vars already in the tracking dir
    for (objName in names(fileMap)) {
        f <- createBindingClosure(objName, trackingEnv)
        makeActiveBinding(objName, env=envir, fun=f)
    }
    setTrackingEnv(trackedEnv=envir, trackingEnv=trackingEnv)
    if (auto) {
        if (opt$debug >= 2)
            cat('track.start: installing track.auto callback\n')
        callbackName <- "track.auto"
        ## remove the old callback (to avoid having duplicates)
        while (is.element(callbackName, getTaskCallbackNames()))
            removeTaskCallback(callbackName)
        addTaskCallback(track.sync.callback, data=envir, name=callbackName)
        assign(".trackAuto", list(on=TRUE, last=-1), envir=trackingEnv)
    }
    if (FALSE && check.Last) {
        ## Stopped using .Last.sys because it was only called when in the globalenv
        if (length(i <- find(".Last.sys")) > 1)
            if (i[1] != find("track.start")[1])
                warning("There are more than one .Last.sys() functions on the search path -- the one from track will is masked and will not run.  This may affect the saving of tracked environments.\n")
            else
                warning("There are more than one .Last.sys() functions on the search path -- the one from track masks others and they will not run\n")
    }
    ## Note that locking the environment is irreversible, and it prevents
    ## rescaning in-place (because the main reason to do that would be to
    ## pick up new variables and delete old ones).  Locking doesn't however
    ## prevent caching, because caching uses the tracking env, not the
    ## tracked env, and the tracking env is not locked.
    if (lockEnv && opt$readonly && environmentName(envir) != "R_GlobalEnv")
        lockEnvironment(envir)
    ## Set up .Last in the globalenv to be track.last(), which will make
    ## sure that all tracked envs are sync'd to disk when R quits.
    ## Do this in the globalenv regardless of which env this call is
    ## tracking, because only .Last in the global env is called when R exits.
    ## This is a good candidate for a different way of doing things, either
    ## a 'Last' hook (doesn't exist in R, but would be nice if it did,
    ## or something like a finalizer on an object, though I wasn't able
    ## to get that to work reliably -- it wasn't always called when R
    ## quitting R.
    .Last <- track.Last
    environment(.Last) <- globalenv()
    existing.Last <- NULL
    ## Fetching an existing .Last here has the desirable side effect that it will be
    ## cached (because .Last is a default member of track.options('alwaysCache'))
    ## Thus, in the case tracking db becomes unavailable, the R-termination
    ## will not be affected by not being able read .Last from disk.
    if (opt$debug >= 2)
        cat('track.start: attempting to install a .Last to save objects to the tracking db\n')
    if (exists(".Last", where=1, inherits=FALSE)) {
        existing.Last <- try(get(".Last", pos=1, inherits=FALSE))
        if (inherits(existing.Last, 'try-error')) {
            if (bindingIsActive(".Last", env=globalenv())) {
                warning(".Last already exists in globalenv as a non-functional active binding -- not installing track.Last; to save objects at end of session, user must call track.stop(all=TRUE) before ending R session (if this copy of .Last is left over from a previous failed initiation of tracking, remove it and try again)")
                existing.Last <- 'broken active.binding'
            } else {
                existing.Last <- 'non-gettable object'
            }
        } else {
            environment(existing.Last) <- globalenv()
        }
    }
    if (!is.null(existing.Last)) {
        if (!identical(.Last, existing.Last) && !identical(existing.Last, 'broken active.binding'))
            warning(".Last already exists in globalenv -- not installing track.Last; to save objects at end of session, user must call track.stop(all=TRUE) before ending R session")
    } else {
        ## There is no existing .Last object in the global environment, so restore the saved one
        Last.pos <- 1
        assign(".Last", .Last, pos=Last.pos)
        ## Do the same thing as in track.sync.callback() for the globalenv
        if (env.is.tracked(pos=1))
            try(track.sync(pos=1, master="envir", taskEnd=TRUE))
    }
    ## Save the tracking summary after working with .Last
    if (!opt$readonly) {
        save.res <- saveObjSummary(trackingEnv, opt=opt, dataDir=getDataDir(dir))
        if (is(save.res, "try-error"))
            stop("could not save '.trackingSummary' in ", attr(save.res, 'file'), ": fix file problem and try again (", save.res, ")")
    }
    assign(".trackingSummaryChanged", FALSE, envir=trackingEnv)
    ## Store the Pid of this R session so that we can identify
    ## situations where a dead .trackEnv has been loaded in by
    ## mistake (probably as a result of saving and reloading an
    ## entire tracked environment.)
    assign(".trackingPid", Sys.getpid(), envir=trackingEnv)
    if (!is.element("track.auto.monitor", getTaskCallbackNames()))
        addTaskCallback(track.auto.monitor, name="track.auto.monitor")
    return(invisible(NULL))
}

track.package.desc <- function(pkg)
    c(paste("Package:",pkg), "Version: 1.0", paste("Date:",date()),
      "Title: Tracked R Objects", "Author: track package", "Maintainer: track package",
      "Description: package of saved objects created by track package", "License: None specified")

# Create a closure that can be used as the active binding
# It has to store the objName and environment where the actual
# object data could be cached.
createBindingClosure <- function(objName, trackingEnv) {
    # Need to force evaluation of the args, otherwise the closure
    # can get the wrong values :-(
    force(objName); force(trackingEnv)
    function(v) {
        if (missing(v))
            getTrackedVar(objName, trackingEnv)
        else
            setTrackedVar(objName, v, trackingEnv)
    }
}
