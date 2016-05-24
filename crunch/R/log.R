getLogBlock <- function (logdf, blockname) {
    ## Extract the subset of the log between BLOCK blockname and the next BLOCK
    blocks <- logdf$scope == "BLOCK"
    blockstart <- blocks & logdf$verb == blockname
    if (any(blockstart)) {
        blockstart <- which(blockstart)[1] + 1## In case there are multiple matches
        blocks <- which(blocks)
        blockend <- blocks > blockstart
        if (any(blockend)) {
            blockend <- blocks[blockend][1] - 1
        } else {
            blockend <- nrow(logdf)
        }
        return(logdf[blockstart:blockend,])
    } else {
        stop("Block ", blockname, " not found", call.=FALSE)
    }
}

blockTimings <- function (logdf) {
    ## Calculate how long each "BLOCK" took
    blocks <- logdf$scope == "BLOCK"
    times <- logdf$timestamp[blocks]
    blocknames <- c("(start)", logdf$verb[blocks])

    return(structure(as.numeric(difftime(c(times, logdf$timestamp[nrow(logdf)]),
        c(logdf$timestamp[1], times), units="secs")),
        .Names=blocknames))
}
