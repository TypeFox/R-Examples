
track.sync.callback <- function(expr, ok, value, visible, data) {
    ## To automatically track new and deleted objects, do
    ##   addTaskCallback(track.sync.callback, data=globalenv())
    ## and
    ##   assign(".trackAuto", list(on=TRUE, last=-1), envir=trackingEnv)
    ## 'data' arg is not used
    trace <- getOption("track.callbacks.trace", FALSE)
    if (trace==1) {
        cat("track.sync.callback: entered at ", date(), ", length(sys.calls())=", length(sys.calls()), "\n", sep="")
        stime <- proc.time()
        on.exit(cat("track.sync.callback: exited at ", date(),
                    " (", paste(round(1000*(proc.time()-stime)[1:3]), c("u", "s", "e"), sep="", collapse=" "), " ms)\n", sep=""))
    }
    ## If we are called from a prompt in a browser, don't do anything
    if (length(sys.calls()) > 1)
        return(TRUE)
    ## Don't repeat the work an explicit call to track.sync()
    ## This caused a crash in R < 2.12 because 'expr' refered to an invalid object
    if (is.call(expr) && as.character(expr[[1]]) == "track.sync")
        return(TRUE)
    ## Loop through all potentially auto-tracked environments and call track.sync where needed
    envs <- search()
    callback.names <- getTaskCallbackNames()
    envs.look <- grep("^(package:|pkgcode:|Autoloads$)", envs, invert=TRUE)
    any.auto.on <- FALSE
    for (i in envs.look) {
        if (env.is.tracked(i) && mget(".trackAuto", ifnotfound=list(list(on=FALSE)),
                                      envir=getTrackingEnv(as.environment(i)))[[1]]$on) {
            any.auto.on <- TRUE
            trackedEnv <- as.environment(i)
            trackingEnv <- getTrackingEnv(trackedEnv, stop.on.not.tracked = FALSE)
            res <- try(track.sync(envir=trackedEnv, trackingEnv=trackingEnv,
                                  full=NA, master="envir", taskEnd=TRUE), silent=TRUE)
            if (is(res, "try-error"))
                warning("oops: track.sync() had a problem (use track.auto(FALSE, pos=) to turn off): ", res)
        }
    }
    if (trace==2) {
        cat("\n")
        flush.console()
    }
    ## Remove the callback when it is no longer wanted
    if (!isTRUE(any.auto.on))
        return(FALSE)
    ## Check that the monitor is alive -- this is a mutual back-scratching exercise
    ## Callbacks can be deleted unintentionally, e.g., by an unfortunately-timed CTRL-C
    if (!is.element("track.auto.monitor", callback.names))
        addTaskCallback(track.auto.monitor, name="track.auto.monitor")
    return(TRUE) # to keep this callback active
}
