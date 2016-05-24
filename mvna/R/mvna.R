mvna <- function(data, state.names, tra, cens.name) {

    if (missing(data))
        stop("Argument 'data' is missing with no default")
    if (missing(tra))
        stop("Argument 'tra' is missing with no default")
    if (missing(state.names))
        stop("Argument 'state.names' is missing with no default")
    if (missing(cens.name))
        stop("Argument 'cens.name' is missing with no default")
    if (!is.data.frame(data))
        stop("Argument 'data' must be a data.frame")
    if (!(xor(sum(c("id", "from", "to", "time") %in% names(data)) != 4,
              sum(c("id", "from", "to", "entry", "exit") %in% names(data)) != 5)))
        stop("'data' must contain the right variables")
    if (nrow(tra) != ncol(tra))
        stop("Argument 'tra' must be a quadratic  matrix.")
    if (sum(diag(tra)) > 0)
        stop("transitions into the same state are not allowed")
    if (nrow(tra) != length(state.names)) {
        stop("The row number of 'tra' must be equal to the number of states.")
    }
    if (!is.logical(tra)) {
        stop("'tra' must be a matrix of logical values, which describes the possible transitions.")
    }
    if (length(state.names) != length(unique(state.names))) {
        stop("The state names must be unique.")
    }
    if (!(is.null(cens.name))) {
        if (cens.name %in% state.names) {
            stop("The name of the censoring variable just is a name of the model states.")
        }
    }
    
### transitions
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
    namen <- paste(trans[, 1], trans[, 2])

    ## test on transitions
    test <- unique(paste(data$from, data$to))
    if (!(is.null(cens.name))) {
        ref <- c(paste(trans$from, trans$to), paste(unique(trans$from), cens.name))
    } else {
        ref <- paste(trans$from, trans$to)
    }
    ref.wo.cens <- paste(trans$from, trans$to)
    if (!(all(test %in% ref)==TRUE))
        stop("There is undefined transitions in the data set")
    if (sum(as.character(data$from)==as.character(data$to)) > 0)
        stop("Transitions into the same state are not allowed")
    if (!(all(ref.wo.cens %in% test) == TRUE))
        warning("You may have specified more possible transitions than actually present in the data")
    
### data.frame transformation
    data$id <- if (is.character(data$id)) as.factor(data$id) else data$id
    data$from <- as.factor(data$from)
    data$to <- as.factor(data$to)
    if (!(is.null(cens.name))) {
        data$from <- factor(data$from, levels = c(cens.name, state.names), ordered = TRUE)
        levels(data$from) <- 0:length(state.names)
        data$to <- factor(data$to, levels = c(cens.name, state.names), ordered = TRUE)
        levels(data$to) <- 0:length(state.names)
    }
    else{
        data$from <- factor(data$from, levels = state.names, ordered = TRUE)
        levels(data$from) <- 1:length(state.names)
        data$to <- factor(data$to, levels = state.names, ordered = TRUE)
        levels(data$to) <- 1:length(state.names)
    }
    
### if not, put like counting process data
    if ("time" %in% names(data)) {
        data <- data[order(data$id, data$time), ]
        idd <- as.integer(factor(data$id))
        entree <- double(length(data$time))
        masque <- rbind(1, apply(as.matrix(idd), 2, diff))
        entree <- c(0, data$time[1:(length(data$time) - 1)]) * (masque == 0)
        data <- data.frame(id = data$id, from = data$from,
                           to = data$to, entry = entree, exit = data$time)
        if (sum(data$entry < data$exit) != nrow(data))
            stop("Exit time from a state must be > entry time")
    } else {
        if (sum(data$entry < data$exit) != nrow(data))
            stop("Exit time from a state must be > entry time")
    }
    
### Computation of the risk set and dN
    ttime <- c(data$entry, data$exit)
    times <- sort(unique(ttime))
    data$from <- as.integer(as.character(data$from))
    data$to <- as.integer(as.character(data$to))
    
    temp <- .C("risk_set_mvna",
               as.integer(nrow(data)),
               as.integer(length(times)),
               as.integer(c(dim(tra), length(times))),
               as.double(times),
               as.integer(data$from),
               as.integer(data$to),
               as.double(data$entry),
               as.double(data$exit),
               nrisk=integer(dim(tra)[1] * length(times)),
               ncens=integer(dim(tra)[1] * length(times)),
               nev=integer(dim(tra)[1] * dim(tra)[2] * length(times)))

    nrisk <- matrix(temp$nrisk, ncol=dim(tra)[1], nrow=length(times))
    ncens <- matrix(temp$ncens, ncol=dim(tra)[1], nrow=length(times))
    nev <- array(temp$nev, dim=c(dim(tra), length(times)))

### computation of the NA estimates
    colnames(nrisk) <- state.names
    dimnames(nev) <- list(state.names, state.names, times)
    
    est <- lapply(1:nrow(trans), function(i) {
        ## on enlève les temps ou nrisk == 0
        ind <- nrisk[, as.character(trans[i, 1])] != 0
        t.nrisk <- nrisk[, as.character(trans[i, 1])][ind]
        t.nev <- nev[as.character(trans[i, 1]), as.character(trans[i, 2]), ][ind]
        na <- cumsum(t.nev/t.nrisk)
        var1 <- cumsum(t.nev/t.nrisk^2)
        var2 <- cumsum(((t.nrisk - t.nev)/t.nrisk^3) * t.nev)
        data.frame(na = na, var.aalen = var1, var.greenwood = var2,
                   time = times[ind])
        })
    
    nrisk <- nrisk[, !(colnames(nrisk) %in% setdiff(unique(trans$to), unique(trans$from))), drop = FALSE]
    names(est) <- namen
    eest <- list(time = times, n.risk = nrisk, n.event = nev, n.cens = ncens,
                 state.names = state.names, cens.name = cens.name,
                 trans = trans)
    res <- c(est, eest)
    
    class(res) <- "mvna"
    res
}    

