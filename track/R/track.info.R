track.info <- function(pos=NULL, envir=as.environment(pos), all=is.null(pos)) {
    if (all) {
        envirs <- search()
        is.tracked <- sapply(envirs, function(envir) env.is.tracked(envir=as.environment(envir)))
        env.list <- lapply(as.list(envirs)[is.tracked], as.environment)
        names(env.list) <- envirs[is.tracked]
        for (j in seq(along=envirs))
            if (!is.tracked[j] && exists(".trackingEnv", envir=as.environment(envirs[j]), inherits=FALSE))
                warning("env ", envirs[j], " (pos ", j, " on search list) appears to be an inactive tracked environment, saved from another session and loaded here inappropriately (see ?track.info)")
        res <- data.frame(row.names=NULL, env.name=envirs[is.tracked],
                          pos=seq(len=length(envirs))[is.tracked])
    } else {
        if (!is.environment(envir))
            envir <- as.environment(envir)
        if (missing(pos))
            if (environmentName(envir)=="R_GlobalEnv")
                pos <- 1
            else
                pos <- match(environmentName(envir), search())
        if (!env.is.tracked(envir=envir))
            stop("env ", envname(envir), " is not tracked")
        res <- data.frame(row.names=NULL, env.name=environmentName(envir), pos=pos)
        env.list <- structure(list(envir), names=environmentName(envir))
    }

    ## This code was for when .Last.sys was a function in the
    ## track package env, but .Last/.Last.sys don't appear
    ## to be run when they live in an attached package.
    if (FALSE && length(i <- find(".Last.sys")) > 1 && !isTRUE(options("global.track.options")$inhibit.Last)) {
        if (i != find("track.start")[1])
            warning("There are more than one .Last.sys() functions on the search path -- the one from track will is masked and will not run.  This may affect the saving of tracked environments.")
        else
            warning("There are more than one .Last.sys() functions on the search path -- the one from track masks others and they will not run.")
    }

    if (nrow(res)) {
        callback.names <- getTaskCallbackNames()
        res <-  cbind(res, data.frame(row.names=NULL, do.call("rbind", lapply(env.list,
                      function(e) {
                           auto <- mget(".trackAuto", ifnotfound=list(list(on=FALSE)),
                                        envir=getTrackingEnv(e))[[1]]$on
                           opt <- track.options(envir=e, c("readonly", "cache"))
                           data.frame(readonly=opt$readonly, cache=opt$cache,
                                      auto=auto, dir=track.dir(envir=e))}
        ))))
        if (any(res$auto) && !is.element("track.auto", callback.names)) {
            pos <- (res$pos[res$auto])[1]
            warning("Task callback track.auto",
                    " is missing; do track.auto(TRUE, pos=", pos, ") to reinstate")
        }
    } else {
        res <- data.frame(env.name=character(0),
                          pos=integer(0),
                          readonly=logical(0),
                          cache=logical(0),
                          auto=logical(0),
                          dir=character(0))
    }
    return(res)
}
