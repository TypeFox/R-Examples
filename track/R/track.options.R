track.options <- function(..., pos=1, envir=as.environment(pos), values=list(...), save=FALSE, clear=FALSE, delete=FALSE, trackingEnv, only.preprocess=FALSE, old.options=list()) {
    ## This function probably tries to do too many things (e.g., in using arg
    ## combinations only.preprocess, trackingEnv, save, old.options, ...)
    ## It would probably be better rewritten into several functions, with
    ##   track.options(..., pos, envir, save) as the user-visible function.
    ##
    ## only.preprocess: return a complete list of new options, with values as specified as arguments,
    ##   and others as defaults
    ## See track.options.Rd DETAILS section for description of valid options
    trackingEnvSupplied <- !missing(trackingEnv) && !is.null(trackingEnv)
    if (only.preprocess) {
        currentOptions <- old.options
    } else {
        if (length(old.options)!=0)
            stop("can only supply old.options when only.preprocess=TRUE")
        if (!trackingEnvSupplied)
            trackingEnv <- getTrackingEnv(envir)
        currentOptions <- mget(".trackingOptions", envir=trackingEnv, ifnotfound=list(list()))[[1]]
        if (length(currentOptions)==0 && exists(".trackingDir", envir=trackingEnv, inherits=FALSE)) {
            ## Read the options from file if we have a trackingDir on envir and no trackingEnv was supplied
            ## I'm not sure if this code ever gets exercised ... the only case I can think of right
            ## now is where .trackingOptions went missing from trackingEnv for some reason.
            ## Try to read them from disk
            dir <- getTrackingDir(envir)
            dataDir <- getDataDir(dir)
            ## we need to work out the RData suffix here -- it's a bit tricky
            ## collect the possibilities and look for them
            gopt <- getOption("global.track.options")
            if (length(gopt$RDataSuffixes)==0)
                gopt$RDataSuffixes <- c("rda", "RData")
            if (!is.character(gopt$RDataSuffixes))
                stop('getOption("global.track.options")$RDataSuffixes must be character data')
            if (any(!regexpr("^[[:alnum:]]+$", gopt$RDataSuffixes)))
                stop('getOption("global.track.options")$RDataSuffixes must consist of alpha-numeric characters only')
            if (length(gopt$RDataSuffixes)==1)
                suffixRegExp <- gopt$RDataSuffixes
            else
                suffixRegExp <- paste("(", paste(gopt$RDataSuffixes, collapse="|", sep=""), ")", sep="")
            suffix <- NULL
            x <- list.files(path=dataDir, pattern=paste("^\\.trackingOptions\\.", suffixRegExp, "$", sep=""), all.files=TRUE)
            if (length(x)>1)
                stop("have multiple options files in '", dataDir, "': ", paste(x, collapse=", "))
            if (length(x)==1) {
                suffix <- sub(".*\\.", "", x)
                file <- file.path(dataDir, paste(".trackingOptions", suffix, sep="."))
                if (!file.exists(file))
                    stop("weird: thought I had the options file, but it doesn't exist...; ", file, "; ", x)
                tmpenv <- new.env(parent=emptyenv())
                load.res <- try(load(file=file, envir=tmpenv), silent=TRUE)
                if (is(load.res, "try-error") || length(load.res)!=1 || load.res!=".trackingOptions") {
                    warning(file, " does not contain a .trackingOptions object -- ignoring it and using system defaults")
                } else {
                    currentOptions <- get(".trackingOptions", envir=tmpenv, inherits=FALSE)
                }
            } else {
                if (trackingEnvSupplied)
                    warning("no .trackingOptions in ", envname(trackingEnv),
                            " and there is no saved options file -- using system defaults")
            }
        }
    }
    if (length(currentOptions)==0) ## in case someone supplied old.options=NULL
        currentOptions <- list()
    ## If we were called like track.options(values=NULL), make this like track.options()
    if (length(values)==1 && is.null(names(values)) && is.null(values[[1]]))
         values <- list()
    optionNames <- c("alwaysCache", "alwaysCacheClass", "alwaysSaveSummary", "autoTrackExcludeClass",
                     "autoTrackExcludePattern", "autoTrackFullSyncWait", "cache",
                     "cacheKeepFun", "cachePolicy", "clobberVars", "compress", "compression_level",
                     "debug", "use.fake.Sys.time", "maintainSummary", "RDataSuffix",
                     "readonly", "stealable", "recordAccesses",
                     "summaryAccess", "summaryTimes", "writeToDisk")
    if (!is.null(names(values))) {
        ## Attempt to set some of the options (including saving to file)
        ## and return the old values.
        ## First retrieve the old values
        query.values <- names(values)
        set.values <- TRUE
        ## No longer a problem:
        ## can't supply trackingEnv and set option values because we wouldn't know where to write the file
        ## if (trackingEnvSupplied)
        ##    stop("cannot supply trackingEnv and set option values")
    } else if (save) {
        query.values <- optionNames
        set.values <- TRUE
    } else {
        ## no names means a query -- expect all char data
        if (length(values)==0) {
            query.values <- optionNames
        } else {
            if (!all(sapply(values, is.character)))
                stop("in a query, all args must be character data")
            query.values <- unlist(values, use.names=FALSE)
        }
        set.values <- FALSE
    }

    if (!all(is.element(query.values, optionNames)))
        stop("unknown option names: ", paste("'", setdiff(query.values, optionNames), "'", sep="", collapse=", "))

    ## See if we need to repair any missing options (shouldn't need to do this)
    ## This is where to set defaults
    if (set.values || only.preprocess)
        need.value <- setdiff(optionNames, names(currentOptions))
    else
        need.value <- setdiff(query.values, names(currentOptions))
    if (length(need.value)) {
        names(need.value) <- need.value
        repaired <- lapply(need.value, function(x)
                           switch(x, cache=TRUE, cachePolicy="eotPurge",
                                  cacheKeepFun="track.plugin.lru",
                                  alwaysCache=c(".Last"),
                                  alwaysCacheClass=c("ff"),
                                  readonly=FALSE,
                                  stealable=FALSE,
                                  writeToDisk=TRUE,
                                  maintainSummary=TRUE, alwaysSaveSummary=FALSE,
                                  recordAccesses=TRUE,
                                  summaryTimes=1, summaryAccess=1, RDataSuffix="rda",
                                  debug=0,
                                  use.fake.Sys.time=FALSE,
                                  autoTrackExcludePattern=c("^\\.track", "^\\.required", "^\\*tmp\\*$", "^.vimplemented", "^.vcoerceable"),
                                  autoTrackExcludeClass=c("RODBC"),
                                  autoTrackFullSyncWait=-1,
                                  clobberVars=c(".Random.seed"),
                                  compress=TRUE, compression_level=1))
        currentOptions <- c(currentOptions, repaired)
    }
    option.values <- currentOptions[query.values]
    if (set.values) {
        new.values <- currentOptions
        for (opt in names(values)) {
            single <- TRUE
            special <- FALSE
            if (opt=="cache") {
                if (!is.logical(values[[opt]]))
                    values[[opt]] <- as.logical(values[[opt]])
            } else if (opt=="cachePolicy") {
                if (!is.character(values[[opt]]))
                    values[[opt]] <- as.character(values[[opt]])
            } else if (opt=="cacheKeepFun") {
                # can be the name of a function, or a function
                # don't do the standard NA & length checks
                special <- TRUE
                f <- values[[opt]]
                if (!identical(f, "none")) {
                    if (!is.function(f)) {
                        if (is.character(f) && length(f)==1) {
                            f <- get(f)
                        } else if (is.name(f)) {
                            f <- eval(f)
                        }
                    }
                    if (!is.function(f) && !is.null(f))
                        stop("cacheKeepFun must be a function or the name of a function")
                    if (!is.null(f) && !all(is.element(c("objs", "inmem", "envname"), names(formals(f)))))
                        stop("cacheKeepFun must have an arguments namd 'objs' and 'envname'")
                }
            } else if (opt=="readonly") {
                if (!is.logical(values[[opt]]))
                    values[[opt]] <- as.logical(values[[opt]])
            } else if (opt=="stealable") {
                if (!is.logical(values[[opt]]))
                    values[[opt]] <- as.logical(values[[opt]])
            } else if (opt=="writeToDisk") {
                if (!is.logical(values[[opt]]))
                    values[[opt]] <- as.logical(values[[opt]])
            } else if (opt=="recordAccesses") {
                if (!is.logical(values[[opt]]))
                    values[[opt]] <- as.logical(values[[opt]])
            } else if (opt=="maintainSummary") {
                if (!is.logical(values[[opt]]))
                    values[[opt]] <- as.logical(values[[opt]])
            } else if (opt=="alwaysSaveSummary") {
                if (!is.logical(values[[opt]]))
                    values[[opt]] <- as.logical(values[[opt]])
            } else if (opt=="summaryTimes") {
                if (!is.integer(values[[opt]]))
                    values[[opt]] <- as.integer(values[[opt]])
            } else if (opt=="summaryAccess") {
                if (!is.integer(values[[opt]]))
                    values[[opt]] <- as.integer(values[[opt]])
            } else if (opt=="RDataSuffix") {
                if (!is.character(values[[opt]]))
                    values[[opt]] <- as.character(values[[opt]])
            } else if (opt=="debug") {
                if (!is.integer(values[[opt]]))
                    values[[opt]] <- as.integer(values[[opt]])
            } else if (opt=="use.fake.Sys.time") {
                if (!is.logical(values[[opt]]))
                    values[[opt]] <- as.logical(values[[opt]])
            } else if (opt=="autoTrackExcludePattern" || opt=="autoTrackExcludeClass"
                       || opt=="alwaysCacheClass") {
                single <- FALSE
                if (!is.character(values[[opt]]))
                    values[[opt]] <- as.character(values[[opt]])
            } else if (opt=="alwaysCache") {
                single <- FALSE
                if (!is.character(values[[opt]]))
                    values[[opt]] <- as.character(values[[opt]])
            } else if (opt=="autoTrackFullSyncWait") {
                if (!is.numeric(values[[opt]]))
                    values[[opt]] <- as.numeric(values[[opt]])
            } else if (opt=="clobberVars") {
                single <- FALSE
                if (!is.character(values[[opt]]))
                    values[[opt]] <- as.character(values[[opt]])
            } else if (opt=="compress") {
                if (!is.character(values[[opt]])) {
                    if (!is.logical(values[[opt]]) || is.na(values[[opt]]))
                        stop("compress must be TRUE/FALSE or the name of a compress technique")
                } else {
                    if (! values[[opt]] %in% c("none", "bzip", "gzip", "xz"))
                        stop("compress as a string must be one of 'none', 'bzip', 'gzip', 'xz'")
                    if (values[[opt]] == 'none')
                        values[[opt]] <- FALSE
                }
            } else if (opt=="compression_level") {
                if (!is.numeric(values[[opt]]) || is.na(values[[opt]]))
                    stop("compression_level must be a number")
                if (values[[opt]] < 1)
                    values[[opt]] <- 1
                if (values[[opt]] > 9)
                    values[[opt]] <- 9
            } else {
                stop("unrecognized option name '", opt, "'")
            }
            if (!special) {
                if (any(is.na(values[[opt]])))
                    stop("cannot set option ", opt, " to an NA value")
                if (single && length(values[[opt]])!=1)
                    stop("option ", opt, " must have a value of length 1")
            }
            ## Now, how we put the value in depends on whether it can have single or multiple values
            ## Do assignment like new.values[opt] <- list(value) so that it works with value==NULL
            if (single || clear)
                new.values[opt] <- list(values[[opt]])
            else if (delete)
                new.values[opt] <- list(setdiff(new.values[[opt]], values[[opt]]))
            else
                new.values[opt] <- list(unique(c(new.values[[opt]], values[[opt]])))
        }

        if (only.preprocess)
            return(new.values)
        if (!new.values$readonly && environmentIsLocked(envir))
            stop("cannot make a readonly tracked environment writable (because cannot unlock a locked environment) -- to make it writeable, use track.detach() followed by track.attach(readonly=FALSE)")
        assign(".trackingOptions", new.values, envir=trackingEnv)
    }
    if (save && !only.preprocess) {
        if (!identical(currentOptions$readonly, FALSE)) {
            warning("cannot save options in a readonly tracking db")
        } else {
            dir <- getTrackingDir(trackingEnv)
            file <- file.path(getDataDir(dir), paste(".trackingOptions", currentOptions$RDataSuffix, sep="."))
            ## if we did change any options, they will have been saved in .trackingOptions in trackingEnv
            save.res <- try(save(list=".trackingOptions", file=file, envir=trackingEnv, compress=FALSE), silent=TRUE)
            if (is(save.res, "try-error"))
                stop("unable to save .trackingOptions in ", file)
        }
    }
    ## Want to return the old values
    if (set.values)
        return(invisible(option.values))
    else
        return(option.values)
}
