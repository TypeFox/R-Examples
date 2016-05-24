### Function to prepare the data in way
### that they can be used in etm()

etmprep <- function(time, status, data, tra, state.names, cens.name = NULL,
                    start = NULL, id = NULL, keep) {

    if (nrow(tra) != ncol(tra))
        stop("'tra' must be quadratic")
    
    ## What are the possible transitions, transient and absorbing states
    if (missing(state.names)) {
        state.names <- as.character(0:(dim(tra)[2] - 1))
    }

    ls <- length(state.names); n <- nrow(data)
    if (ls != dim(tra)[2])
        stop("Discrepancy between 'tra' and the number of states specified in 'state.names'")

    if (length(time) != ls) {
        stop("The length of 'time' must be equal to the number of states")
    }
    
    colnames(tra) <- rownames(tra) <- state.names
    t.from <- lapply(1:dim(tra)[2], function(i) {
        rep(rownames(tra)[i], sum(tra[i, ]))
    })
    t.from <- unlist(t.from)
    t.to <- lapply(1:dim(tra)[2], function(i) {
        colnames(tra)[tra[i, ]==TRUE]
    })
    t.to <- unlist(t.to)
    trans <- data.frame(from=t.from, to=t.to)
    absorb <- setdiff(levels(trans$to), levels(trans$from))
    transient <- unique(state.names[!(state.names %in% absorb)])

    ## extract informations in time
    ind <- match(time[!is.na(time)], names(data))
    if (any(is.na(ind)))
        stop("At least one element in 'time' is not in 'data'")
    indd <- which(time %in% names(data))
    time <- matrix(NA, n, ls)
    time[, indd] <- as.matrix(data[, ind])

    ## extract infos in status
    if (length(status) != ls) {
        stop("The length of 'status' must be equal to the number of states")
    }
    ind <- match(status[!is.na(status)], names(data))
    if (any(is.na(ind)))
        stop("At least one element in 'status' is not in 'data'")
    indd <- which(status %in% names(data))
    status <- matrix(NA, n, ls)
    status[, indd] <- as.matrix(data[, ind])

    if (is.null(start)) {
        start.state <- rep(state.names[1], n)
        start.time <- rep(0, n)
    } else {
        if ((length(start$state) != nrow(data)) | (length(start$time) != nrow(data))) 
            stop("'start$state' or 'start$time' are not as long as the data")
        if (!all(unique(start$state) %in% state.names))
            stop("'start$state' not in 'state.names'")
        start.state <- start$state
        start.time <- start$time
    }

    if (is.null(id)) {
        id <- seq_len(n)
    } else id <- data[, id]

    if (!missing(keep)) {
        cova <- data[, keep, drop = FALSE]
    } else keep <- NULL

    ## let's try to start the real work
    newdata <- lapply(seq_len(n), function(i) {
        ind <- which(status[i, ] != 0)
        li <- length(ind)
        if (li == 0) {
            from <- start.state[i]
            to <- cens.name
            entry <- start.time[i]
            exit <- time[i, ncol(time)]
            idd <- id[i]
        } else {
            from <- c(start.state[i], state.names[ind[-li]])
            to <- state.names[ind]
            entry <- c(start.time[i], time[i, ind[-li]])
            exit <- time[i, ind]
            idd <- rep(id[i], length(exit))
            if (to[length(to)] %in% transient) {
                from <- c(from, to[length(to)])
                to <- c(to, cens.name)
                entry <- c(entry, exit[length(exit)])
                exit <- c(exit, time[i, ncol(time)])
                idd <- c(idd, id[i])
            }
        }

        if (is.null(keep)) {
            tmp <- data.frame(idd, entry, exit, from, to)
        } else {
            aa <- matrix(apply(cova[i, , drop = FALSE], 2, rep, length(exit)),
                         length(exit), ncol(cova))
            tmp <- data.frame(idd, entry, exit, from, to, aa)
        }
        tmp
    })
    newdata <- do.call(rbind, newdata)
    names(newdata) <- c("id", "entry", "exit", "from", "to", keep)
    if (is.factor(newdata$from) || is.factor(newdata$to)) {
        aa <- unique(c(levels(newdata$from), levels(newdata$to)))
        newdata$from <- factor(as.character(newdata$from), levels = aa)
        newdata$to <- factor(as.character(newdata$to), levels = aa)
    }
   
    newdata
}
