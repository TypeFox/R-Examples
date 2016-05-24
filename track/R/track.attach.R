track.attach <- function(dir, pos=2, name=NULL, create=FALSE, readonly=!create, lockEnv=FALSE, verbose=TRUE, auto=NULL, dup.ok=FALSE) {
    ## if (missing(dir) && missing(pos) && is.null(name) && !readonly && !lockEnv)
    ##    return(track.start())
    if (pos < 2)
        stop("pos must be >= 2")
    if (file.exists(file.path(dir, "filemap.txt"))) {
        trackingDir <- dir
        if (basename(dir)=="rdatadir")
            dir <- dirname(dir)
    } else if (file.exists(file.path(dir, "rdatadir", "filemap.txt"))) {
        trackingDir <- file.path(dir, "rdatadir")
    } else {
        if (!create)
            stop("tracking db does not exist in '", dir, "' or '", file.path(dir, "rdatadir"),
                 "' and cannot create because create=FALSE")
        else if (readonly)
            stop("tracking db does not exist in '", dir, "' or '", file.path(dir, "rdatadir"),
                 "' and cannot create because readonly=TRUE")
        trackingDir <- dir
    }
    ## dir is the parent of trackingDir
    abs.dir <- getAbsolutePath(dir)
    dup.count <- 0
    if (dir.exists(abs.dir)) {
        ## normalize.path only works if the file exists
        rel.dir <- quietNormalizePath(find.relative.path(getwd(), abs.dir), mustWork=FALSE, winslash='/')
        search.list <- search()
        if (!is.na(rel.dir)) for (i in seq(along=search.list)) {
            i.env <- as.environment(i)
            if (env.is.tracked(envir=i.env)) {
                i.path <- quietNormalizePath(find.relative.path(getwd(), track.datadir(envir=i.env)),
                                        mustWork=FALSE, winslash='/')
                if (!is.na(i.path) && rel.dir==i.path) {
                    if (dup.ok) {
                        dup.count <- dup.count + 1
                        if (dup.count==1 && dup.ok < 2)
                            warning('reattaching already attached tracking database ', rel.dir)
                    } else {
                        warning('tracking database ', rel.dir,
                                ' is already attached at position ', i,
                                ' on search list; not reattaching')
                        return(invisible(NULL))
                    }
                }
            }
        }
    }
    if (is.null(name))
        if (dup.count > 0)
            name <- paste(abs.dir, '[', dup.count+1, ']', sep='')
        else
            name <- abs.dir
    attach(what=NULL, pos=pos, name=name)
    assign(".trackingCreated", TRUE, pos=pos)
    if (verbose)
        cat("Attaching tracking db in '", dir, "' to env in pos ", pos,
            if (readonly) " (readonly)" else " (writable)", "\n", sep="")
    return(track.start(trackingDir, pos=pos, readonly=readonly,
                       create=create, lockEnv=lockEnv, check.Last=FALSE,
                       verbose=FALSE, auto=auto))
}
