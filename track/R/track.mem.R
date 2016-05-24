track.mem <- function(pos=NULL, envir=as.environment(pos), all=is.null(pos)) {
    if (all) {
        envirs <- search()
        is.tracked <- sapply(envirs, function(envir) env.is.tracked(envir=as.environment(envir)))
        env.list <- lapply(as.list(envirs)[is.tracked], as.environment)
        names(env.list) <- envirs[is.tracked]
        for (j in seq(along=envirs))
            if (!is.tracked[j] && exists(".trackingEnv", envir=as.environment(envirs[j]), inherits=FALSE))
                warning("env ", envirs[j], " (pos ", j, " on search list) appears to be an inactive tracked environment, saved from another session and loaded here inappropriately (see ?track.info)")
        res <- do.call('rbind', lapply(seq(len=length(envirs))[is.tracked], function(pos) track.mem(pos=pos, all=FALSE)))
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
        trackingEnv <- getTrackingEnv(envir)
        # objs <- track.summary(envir=envir, all=TRUE)
        objs <- get('.trackingSummary', envir=trackingEnv)
        objs$InMem <- is.element(rownames(objs), ls(envir=trackingEnv, all.names=TRUE))
        res <- data.frame(pos=pos, envName=search()[pos], nInMem=sum(objs$InMem, na.rm=TRUE),
                          MBInMem=round(sum(objs$size[objs$InMem], na.rm=TRUE)/1e6, 3),
                          nOutMem=sum(!objs$InMem, na.rm=TRUE),
                          MBOutMem=round(sum(objs$size[!objs$InMem], na.rm=TRUE)/1e6, 3))
    }
    return(res)
}
