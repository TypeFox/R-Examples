track.auto.monitor <- function(expr, value, ok, visible) {
    ## A callback that monitors existence of task callbacks
    ## for track.sync and re-adds them if they have gone
    ## missing (task callbacks can go missing if they stop
    ## with an error, which is what happens when the user
    ## interrupts the R process, or if there are transient
    ## file system problems.)  Try to be fast -- only look
    ## at environments on the search path that are
    ## candidates for tracking (i.e., not packages).
    trace <- getOption("track.callbacks.trace", FALSE)
    if (trace==1) {
        cat("track.auto.monitor: entered at ", date(), "\n", sep="")
        stime <- proc.time()
    }
    if (identical(getOption("incr.hist.active", FALSE), TRUE)) {
        if (!is.element("track.history.writer", getTaskCallbackNames())) {
            cat("track.auto.monitor: Task callback track.history.writer seems to have disappeared; reinstating...\n")
            addTaskCallback(track.history.writer, name="track.history.writer")
        }
    }
    envs <- search()
    callback.names <- getTaskCallbackNames()
    envs.look <- grep("^(package:|pkgcode:|Autoloads$)", envs, invert=TRUE)
    any.auto.on <- FALSE
    for (i in envs.look) {
        if (env.is.tracked(i)) {
            if (mget(".trackAuto", ifnotfound=list(list(on=FALSE)),
                     envir=getTrackingEnv(as.environment(i)))[[1]]$on) {
                any.auto.on <- TRUE
                pos <- i
            }
        }
    }
    if (any.auto.on) {
        callback.name <- "track.auto"
        if (!is.element(callback.name, callback.names)) {
            cat("track.auto.monitor: Task callback", callback.name,
                "seems to have disappeared; reinstating...\n")
            try(track.auto(TRUE, pos = pos))
        }
    }
    if (trace==1) {
        cat("track.auto.monitor: exited at ", date(),
             " (", paste(round(1000*(proc.time()-stime)[1:3]),
                         c("u", "s", "e"), sep="", collapse=" "), ")\n", sep="")
    }
    TRUE
}
