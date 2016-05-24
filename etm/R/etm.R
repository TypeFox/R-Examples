prodint <- function(dna, times, first, last, indi) {
    I <- array(0, dim=dim(dna)[c(1, 2)])
    diag(I) <- 1
    if (first >= last) {
        est <- array(I, dim=c(dim(dna)[c(1, 2)], 1))
        time <- NULL
    } else {
        est <- array(0, dim=c(dim(dna)[c(1, 2)], (last-first+1)))
        est[, , 1] <- I + dna[, , first] * indi[1]
        j <- 2
        for (i in (first + 1):last) {
            est[, , j] <- est[, , j-1] %*% (I + dna[, , i] * indi[j])
            j <- j + 1
        }
        time <- times[first:last]
    }
    list(est=est, time=time)
}



#################################################
### Variance Lai and Ying for competing risks ###
#################################################

var.ly <- function(est, state.names, nrisk, nev, times, first, last, indi) {

    if (first >= last) {
        return(NULL)
        
    } else {

        nCompRisks <- length(state.names) - 1
                
        ## prepare what we need
        cif <- n.event <- matrix(nrow = last - first + 1, ncol = nCompRisks)
        nev <- nev[, , first:last]
        nrisk <- nrisk[first:last, ]
        time <- times[first:last]
        lt <- length(time)

        for (i in seq_len(last - first + 1)) {
            cif[i, ] <- est[1, 2:(nCompRisks + 1), i]
            n.event[i, ] <- nev[1, 2:(nCompRisks + 1), i]
        }
        sminus <- c(1, est[1, 1, 1:(last - first)])
        S <- est[1, 1, ]

        ## create the matrix of covariances
        out <- array(0, dim = c((nCompRisks + 1)^2, (nCompRisks + 1)^2, (last-first+1)))

        ## get the indices on where to put the variance
        pos <- sapply(1:length(state.names), function(i) {
            paste(state.names, state.names[i])
        })
        pos <- matrix(pos)
        dimnames(out) <- list(pos, pos, time)
        pos.cp <- sapply(seq_along(state.names), function(i)
            paste(state.names[1], state.names[i], sep = " "))[-1]
        ind.cp <- which(pos %in%  pos.cp, arr.ind = TRUE)

        ## the real shebang
        for (i in seq_along(ind.cp)) {
            for (j in seq_len(lt)) {
                f <- cif[1:j, i]
                s <- sminus[1:j]
                spasminus <- S[1:j]
                y <- nrisk[1:j, 1]
                dn <- rowSums(array(n.event[1:j, ], dim = c(j, nCompRisks)))
                dnt <- n.event[1:j, i]
                indi.loop <- indi[1:j]
                ## from biomJ paper eq. (6)
                vly <- sum(((f[j] - f)^2 / (y - dn)) * (dn/y) * indi.loop +
                    s^2/y^3 * (y - dnt - 2 * (y - dn) * ((f[j] - f)/spasminus)) * dnt * indi.loop)
                out[ind.cp[i], ind.cp[i], j] <- vly
            }
        }
    }

    return(out)
}


####################################
### Variance of the AJ estimator ###
####################################

var.aj <- function(est, dna, nrisk, nev, times, first, last) {
    d <- dim(nev)[1]
    if (first >= last) {
        return(NULL)
    } else {
        out <- array(0, dim=c(dim(dna)[c(1, 2)]^2, (last-first+1)))
        cov.dna <- matrix(.C(cov_dna,
                             as.double(nrisk[first, ]),
                             as.double(nev[, , first]),
                             as.integer(d),
                             cov = double(d^2 * d^2)
                             )$cov, d^2, d^2)
        bI <- diag(1, d^2)
        out[, , 1] <- bI %*% cov.dna %*% bI
        Id <- diag(1, d)
        for (i in 1:length(times[(first + 1):last])) {
            step <- first + i
            cov.dna <- matrix(.C(cov_dna,
                                 as.double(nrisk[step, ]),
                                 as.double(nev[, , step]),
                                 as.integer(d),
                                 cov = double(d^2 * d^2)
                                 )$cov, d^2, d^2)
            out[, , i + 1] <- (t(Id + dna[, , step]) %x% Id) %*% out[, , i] %*%
                ((Id + dna[, , step]) %x% Id) +
                    (Id %x% est[, , i]) %*% cov.dna %*% (Id %x% t(est[, , i]))
        }
    }
    return(out)
}
         


###########
### etm ###
###########      
      
etm <- function(data, state.names, tra, cens.name, s, t="last",
                covariance=TRUE, delta.na = TRUE, modif = FALSE,
                alpha = 1/4, c = 1) {
    
    if (missing(data))
        stop("Argument 'data' is missing with no default")
    if (missing(tra))
        stop("Argument 'tra' is missing with no default")
    if (missing(state.names))
        stop("Argument 'state.names' is missing with no default")
    if (missing(cens.name))
        stop("Argument 'cens.name' is missing with no default")
    if (missing(s))
        stop("Argument 's' is missing with no default")
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

    ## if modif TRUE, check that the model is competing risks. else
    ## set to false and issue a warning
    if (modif == TRUE && covariance == TRUE) {
        ## check for competing risks
        tr.cp <- tra_comp(length(state.names) - 1)
        if (any(dim(tra) != dim(tr.cp)) | (all(dim(tra) == dim(tr.cp)) && !all(tra == tr.cp))) {
            covariance <- FALSE
            warning("The variance of the estimator with the Lay and Ying transformation is only computed for competing risks data")
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

    n <- length(unique(data$id))
### data.frame transformation
    data$id <- if (is.character(data$id)) as.factor(data$id) else data$id
    data$from <- as.factor(data$from)
    data$to <- as.factor(data$to)
    if (!(is.null(cens.name))) {
        data$from <- factor(data$from, levels = c(cens.name, state.names), ordered = TRUE)
        levels(data$from) <- 0:length(state.names)
        data$to <- factor(data$to, levels = c(cens.name, state.names), ordered = TRUE)
        levels(data$to) <- 0:length(state.names)
    } else{
        data$from <- factor(data$from, levels = state.names, ordered = TRUE)
        levels(data$from) <- 1:length(state.names)
        data$to <- factor(data$to, levels = state.names, ordered = TRUE)
        levels(data$to) <- 1:length(state.names)
    }
    
### if not, put like counting process data
    if ("time" %in% names(data)) {
        data <- data[order(data$id, data$time), ]
        idd <- as.integer(data$id)
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
    temp <- .C(risk_set_etm,
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
               nev=integer(dim(tra)[1] * dim(tra)[2] * length(times)),
               dna=double(dim(tra)[1] * dim(tra)[2] * length(times)))
    
    nrisk <- matrix(temp$nrisk, ncol=dim(tra)[1], nrow=length(times))
    ncens <- matrix(temp$ncens, ncol=dim(tra)[1], nrow=length(times))
    nev <- array(temp$nev, dim=c(dim(tra), length(times)))
    dna <- array(temp$dna, dim=c(dim(tra), length(times)))
    ii <- seq_len(dim(tra)[1])
    for (i in seq_along(times)) {
        dna[cbind(ii, ii, i)] <- -(.rowSums(nev[, , i], dim(nev)[1], dim(nev)[1], FALSE))/nrisk[i, ]
        ## dna[cbind(ii, ii, i)] <- -(rowSums(nev[, , i])/nrisk[i, ])
    }
    dna[is.nan(dna)] <- 0
    
### computation of the Aalen-Johansen estimator
    if (t=="last") t <- times[length(times)]
    if (!(0 <= s & s < t))
        stop("'s' and 't' must be positive, and s < t")
    if (t <= times[1] | s >= times[length(times)])
        stop("'s' or 't' is an invalid time")
    first <- length(times[times <= s]) + 1
    last <- length(times[times <= t])
    
    if (first >= last) {
        est <- list()
        est$est <- array(diag(1, dim(tra)[1], dim(tra)[2]), c(dim(tra), 1))
        dimnames(est$est) <- list(state.names, state.names, t)
        est$time <- NULL
        var <- NULL
        nrisk <- matrix(nrisk[last, ], 1, dim(tra)[1])
        nev <- array(0, dim(tra))
        
    } else {
        
        aa <- nrisk[first:last, ]
        if (modif) {
            which.compute <- as.integer(aa >= c * n^alpha)
        } else {
            which.compute <- rep(1, length(aa))
        }
        est <- prodint(dna, times, first, last, which.compute)
        
        if (covariance == TRUE) {
            if (modif == FALSE) {
                var <- var.aj(est$est, dna, nrisk, nev, times, first, last)
                pos <- sapply(1:length(state.names), function(i) {
                    paste(state.names, state.names[i])
                })
                pos <- matrix(pos)
                dimnames(var) <- list(pos, pos, est$time)
                
            }  else {
                
                var <- var.ly(est$est, state.names, nrisk, nev, times, first, last, which.compute)
            }
            
        } else {
            var <- NULL
        }
        
        if (delta.na) {
            delta.na <- dna[, , first:last]
        }
        else delta.na <- NULL
        
        nrisk <- nrisk[first:last, ]
        nev <- nev[, , first:last]
        dimnames(est$est) <- list(state.names, state.names, est$time)
        dimnames(nev) <- list(state.names, state.names, est$time)
    }
    
    colnames(nrisk) <- state.names
    nrisk <- nrisk[, !(colnames(nrisk) %in%
                       setdiff(unique(trans$to), unique(trans$from))),
                   drop = FALSE]
    
    res <- list(est = est$est, cov = var, time = est$time, s =s, t = t,
                trans = trans, state.names = state.names,
                cens.name = cens.name,
                n.risk = nrisk, n.event = nev, delta.na = delta.na,
                ind.n.risk = ceiling(c * n^alpha))
    class(res) <- "etm"
    res
}

