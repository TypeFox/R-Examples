trprob <- function(x, ...) {
    UseMethod("trprob")
}

trcov <- function(x, ...) {
    UseMethod("trcov")
}

trprob.etm <- function(x, tr.choice, timepoints, ...) {
    if (!inherits(x, "etm"))
        stop("'x' must be a 'etm' object")
    if (!is.character(tr.choice))
        stop("'tr.choice' must be a character vector")
    if (length(tr.choice) != 1)
        stop("The function only extracts 1 transition probability")
    pos <- sapply(1:length(x$state.names), function(i) {
        paste(x$state.names, x$state.names[i])
    })
    pos <- matrix(pos)
    if (!(tr.choice %in% pos))
        stop("'tr.choice' not in the possible transitions")
    trans.sep <- strsplit(tr.choice, " ")
    if (length(trans.sep[[1]]) != 2) {
        tt <- charmatch(trans.sep[[1]], x$state.names, nomatch = 0)
        trans.sep[[1]] <- x$state.names[tt]
    }
    trans.sep <- unlist(trans.sep)

    if (missing(timepoints)) {
        tmp <- x$est[trans.sep[1], trans.sep[2], ]
    }
    else {
        ind <- findInterval(timepoints, x$time)
        tmp <- numeric(length(timepoints))
        place <- which(ind != 0)
        tmp[place] <- x$est[trans.sep[1], trans.sep[2], ind]
    }
    tmp
}
    
trcov.etm <- function(x, tr.choice, timepoints, ...) {
    if (!inherits(x, "etm"))
        stop("'x' must be a 'etm' object")
    if (!is.character(tr.choice))
        stop("'tr.choice' must be a character vector")
    if (!(length(tr.choice) %in% c(1, 2)))
        stop("'tr.choice' must be of length 1 or 2")
        pos <- sapply(1:length(x$state.names), function(i) {
        paste(x$state.names, x$state.names[i])
    })
    pos <- matrix(pos)
    if (!all((tr.choice %in% pos)))
        stop("'tr.choice' not in the possible transitions")
    if (length(tr.choice) == 1) {
        tr.choice <- rep(tr.choice, 2)
    }
    if (missing(timepoints)) {
        tmp <- x$cov[tr.choice[1], tr.choice[2], ]
    }
    else {
        ind <- findInterval(timepoints, x$time)
        tmp <- numeric(length(timepoints))
        place <- which(ind != 0)
        tmp[place] <- x$cov[tr.choice[1], tr.choice[2], ind]
    }
    tmp
}
