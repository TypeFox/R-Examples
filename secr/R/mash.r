# 2012-09-17 moved from trap.builder

mash <- function(object, origin = c(0,0), clustergroup = NULL, ...) {

## mash() recasts a capthist object in which the detectors belong to
## multiple clusters as a capthist with multiple detections at one cluster
## This assumes independence of clusters: if any individuals were detected
## on multiple clusters their new detection histories will be misleading

    if (is.list(clustergroup) & (length(clustergroup) > 1)) {
        if (ms(object))
            stop ("cannot regroup multisession capthist")
        out <- vector('list')
        for (i in 1:length(clustergroup)) {
            out[[i]] <- mash (object, origin, clustergroup[[i]], ...)
        }
        names(out) <- names(clustergroup)
        class(out) <- c('list', 'capthist')
        if (length(out) == 1) out <- out[[1]]
        return(out)
    }
    else if (ms(object)) {
        out <- lapply(object, mash, origin, clustergroup, ...)
#        names(out) <- names(clustergroup)
        names(out) <- names(object)
        class(out) <- c('list', 'capthist')
        if (length(out) == 1) out <- out[[1]]
        return(out)
    }
    else {
        if (!is.null(clustergroup)) {
            trapsi <- clusterID(traps(object)) %in% clustergroup
            object <- subset(object, traps = trapsi)
        }
        trps <- traps(object)
        if (!is.null(covariates(trps)))
            warning ("detector covariates are discarded by mash()")
        if (!is.null(usage(trps)))
            warning ("usage discarded by mash()")
        cluster <- clusterID(trps)
        centres <- cluster.centres(trps)

        ## how many individuals per cluster?
        ## assign each to the first cluster in which it appears
        cl <- cluster[trap(object, names = FALSE)]
        ID <- animalID(object, names = FALSE)
        n.mash <- table (cl[match(unique(ID),ID)])
        class(n.mash) <- 'integer'  ## from 'table'

        if (is.null(cluster))
            stop ("requires cluster covariate")
        tmp <- split(trps, cluster)
        if (length(unique(sapply(tmp, nrow))) != 1)
            warning ("unequal number of detectors per cluster")

        trapnum <- clustertrap(trps)
        if (is.null(trapnum)) {
            tmp <- lapply(tmp, function(x) {x$trapnum <- 1:nrow(x); x})
            trapnum <- unlist(sapply(tmp, function(x) x$trapnum))
        }

        ## take first cluster for new traps
        newtraps <- tmp[[1]]
        rownames(newtraps) <- 1:nrow(newtraps)
        mxy <- apply(newtraps, 2, min)
        newtraps <- shift(newtraps, origin-mxy[1:2])

        ## added 2012=07-26
        attr(newtraps, 'cluster') <- NULL
        attr(newtraps, 'clustertrap') <- NULL
        attr(newtraps, 'covariates') <- NULL

        sigcov <- NULL
        if ( length(animalID(object)) == 0) {
            tempdf <- data.frame(
                session = session(object),
                ID = 'NONE',
                occ = ncol(object),
                trap = 1)
        }
        else {
            tempdf <- data.frame(
                session = rep(session(object), length(animalID(object))),
                ID = animalID(object),
                occ = occasion(object),
                trap = trapnum[trap(object, names=FALSE)]
            )
            if (!is.null(attr(object, 'signalframe'))) {
                tempdf <- cbind(tempdf, attr(object, 'signalframe'))
                sigcov <- names(tempdf)[!(names(tempdf) %in% c('signal','noise'))]
            }
        }
        tempcapt <- make.capthist(tempdf, newtraps, cutval = attr(object, "cutval"),
                                  signalcovariates = sigcov, ...)

        attr(tempcapt, 'n.mash') <- n.mash
        attr(tempcapt, 'centres') <- centres
        tempcapt
    }
}

