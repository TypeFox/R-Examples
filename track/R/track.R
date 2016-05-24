track <- function(expr, pos=1, envir=as.environment(pos), list=NULL, pattern=NULL, glob=NULL, exclude=TRUE) {
    trackingEnv <- getTrackingEnv(envir)
    opt <- track.options(trackingEnv=trackingEnv)
    if (opt$readonly)
        stop(envname(trackingEnv), " is readonly")
    haveVal <- FALSE
    if (!missing(expr)) {
        ## evaluate expr if necessary, and convert to list
        qexpr <- substitute(expr)
        if (is.name(qexpr) || is.character(qexpr)) {
            objName <- as.character(qexpr)
        } else if (mode(qexpr)=="call" && class(qexpr)=="<-") {
            if (length(list))
                stop("cannot use both assignment expr and list at the same time in track(LHS <- RHS, list=...)")
            if (!is.name(qexpr[[2]]))
                stop("LHS must be a simple var name in track(LHS <- RHS)")
            ## need to evaluate and assign
            objName <- as.character(qexpr[[2]])
            objval <- eval(qexpr[[3]], envir=parent.frame(), enclos=parent.frame())
            haveVal <- TRUE
            ## check if the var already exists
            if (exists(objName, envir=envir, inherits=FALSE)) {
                if (objIsTracked(objName, envir, trackingEnv)) {
                    setTrackedVar(objName, objval, trackingEnv, opt)
                    ## will detect that this var is already tracked below, and won't
                    ## do setTrackedVar() again
                } else {
                    remove(list=objName, envir=envir)
                }
            }
        } else {
            stop("argument to track() must be a quoted or unquoted variable name, or an assignment")
        }
        list <- c(objName, list)
    } else if (is.null(list)) {
        list <- untracked(envir=envir, pattern=pattern, glob=glob)
        if (isTRUE(exclude))
            exclude <- opt$autoTrackExcludePattern
        if (identical(exclude, FALSE))
            exclude <- NULL
        for (re in exclude)
            list <- grep(re, list, invert=TRUE, value=TRUE)
    }
    if (length(list)) {
        alreadyTracked <- objIsTracked(list, envir, trackingEnv)
        if (any(alreadyTracked)) {
            # this is a pretty useless warning, so skip it
            # warning("the following objects are already tracked: ",
            #          paste("'", list[which(alreadyTracked)[seq(len=min(3,sum(alreadyTracked)))]], "'", sep="", collapse=", "),
            #                if (sum(alreadyTracked) > 3) ", ...")
            list <- list[!alreadyTracked]
        }
    }
    i <- which(isReservedName(list))
    if (length(i)) {
        warning("cannot track ", length(i), " variables (these are used in implementing tracking or are illegal variable names): ",
                paste("'", list[i], "'", sep="", collapse=", "))
        list <- list[-i]
    }
    all.objs <- ls(envir=envir, all.names=TRUE)
    for (objName in list) {
        ## Doesn't matter if it already exists....
        ## if (exists(objName, trackingEnv))
        ##     stop("'", objName, "' already exists in trackingEnv ", envname(trackingEnv))
        if (haveVal) {
            ## Nothing to do here -- already have the val in objval
        } else if (!is.element(objName, all.objs)) {
            objval <- NULL
        } else {
            if (bindingIsActive(objName, envir)) {
                warning("cannot track '", objName, "' because it is an active binding")
                next
            }
            objval <- get(objName, envir=envir, inherits=FALSE)
            remove(list=objName, envir=envir)
            ## robustness danger point: can have objval in here, but
            ## assigned in no environment -- might be better to remove
            ## it only after setTrackedVar() has succeeded.  However,
            ## setTrackedVar cannot work if it is still here...
        }
        ## robustness: what to do if the assign inside setTrackedVar fails?
        setTrackedVar(objName, objval, trackingEnv, opt)
        f <- createBindingClosure(objName, trackingEnv)
        makeActiveBinding(objName, env=envir, fun=f)
    }
    return(invisible(list))
}
