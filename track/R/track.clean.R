track.clean <- function() {
    ## Not sure that I want this function; hence it is not exported.
    ## Loop through all tracked environments and call track.sync
    ## with full=TRUE, which will do a potentially time-consuming
    ## check that all tracked variables have an activeBinding.
    ## It is possible for a tracked object to loose its activeBinding
    ## when a variable is removed and reassigned within one top
    ## level call (merely replacing it would go via the activeBinding,
    ## with out messing anything up.)
    envs <- search()
    callback.names <- getTaskCallbackNames()
    envs.look <- grep("^(package:|pkgcode:|Autoloads$)", envs, invert=TRUE)
    for (i in envs.look) {
        if (env.is.tracked(i)) {
            trackedEnv <- as.environment(i)
            trackingEnv <- getTrackingEnv(trackedEnv, stop.on.not.tracked = FALSE)
            res <- try(track.sync(envir=trackedEnv, trackingEnv=trackingEnv,
                                  full=TRUE, master="envir", taskEnd=TRUE), silent=TRUE)
            if (is(res, "try-error"))
                warning("oops: track.sync() had a problem: ", res)
        }
    }
}
