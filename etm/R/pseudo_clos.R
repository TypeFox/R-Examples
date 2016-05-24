### Function to compute the pseudo values
## Modelling will be done in another function to offer more
## flexibility
### Author: Arthur Allignol <arthur.allignol@fdm.uni-freiburg.de>

closPseudo <- function(data, state.names, tra, cens.name, s = 0,
                       formula, aw = FALSE, ratio = FALSE, ncores = 1) {

    ## take care of the formula argument
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "id", "subset", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    Terms <- if (missing(data)) terms(formula)
    else terms(formula, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())

    ids <- unique(data$id)
    n <- length(ids)

    ## theta. From there we'll see what kind of model it is
    ## is no alternative weights, NULL
    ## No competing risks: not in the list
    theta <- unlist(clos(etm(data = data, state.names = state.names, tra = tra,
                             cens.name = cens.name, s = 0, covariance = FALSE),
                         aw = aw, ratio = ratio)[c("e.phi", "e.phi.weights.1",
                         "e.phi.weights.other",
                         "e.phi2", "e.phi3")])

    competing <- "e.phi2" %in% names(theta)

    ## Compute pseudo values, and store results depending of competing
    ## and aw
    ## TODO: ACTUALLY COMPUTE THE PSEUDO VALUES
    namen <- c("ps.e.phi", "ps.e.phi.weights.1", "ps.e.phi.weights.other",
               "ps.e.phi2", "ps.e.phi3")

    psMatrix <- parallel::mclapply(seq_along(ids), function(i) {
        temp <- clos(etm(data = data[!(data$id %in% ids[i]), ],
                         state.names = state.names, tra = tra,
                         cens.name = cens.name, s = 0, covariance = FALSE),
                     aw = aw, ratio = ratio)
        
        cbind(temp$e.phi, temp$e.phi.weights.1, temp$e.phi.weights.other,
              temp$e.phi2, temp$e.phi3)
    }, mc.cores = ncores)
    ##}  else {
    ##       psMatrix <- lapply(seq_along(ids), function(i) {
    ##           temp <- clos(etm(data = data[!(data$id %in% ids[i]), ],
    ##                            state.names = state.names, tra = tra,
    ##                            cens.name = cens.name, s = 0, covariance = FALSE),
    ##                        aw = aw, ratio = ratio)
    
    ##           cbind(temp$e.phi, temp$e.phi.weights.1, temp$e.phi.weights.other,
    ##                 temp$e.phi2, temp$e.phi3)
    ##       })
    ##   }
             
    psMatrix <- data.frame(do.call(rbind, psMatrix))
    
    psMatrix <- lapply(seq_along(psMatrix), function(i) {
        n * theta[i] - (n - 1) * psMatrix[, i]
    })
    psMatrix <- do.call(cbind, psMatrix)
    colnames(psMatrix) <- namen[c(TRUE, aw, aw, competing, competing)]
    ## the pseudo values n * ref - (n - 1) * temp
    ## psMatrix <- matrix(apply(psMatrix, 1, function(x) n * theta - (n - 1) * x),
    ##                    nrow = dim(psMatrix)[1], ncol = dim(psMatrix)[2])
    ## colnames(psMatrix) <- namen[c(TRUE, aw, aw, competing, competing)]
    
    cov <- m[!duplicated(data$id), , drop = FALSE]
    colnames(cov) <- attr(Terms, "term.labels")
    
    theta <- matrix(theta, nrow = 1)
    colnames(theta) <- c("e.phi", "e.phi.weights.1",
                         "e.phi.weights.other", "e.phi2",
                         "e.phi3")[c(TRUE, aw, aw, competing, competing)]
        
    zzz <- list(pseudoData = data.frame(id = ids, psMatrix, cov),
                theta = theta, aw = aw, call = call)
    class(zzz) <- "closPseudo"

    zzz
}


### A function to compute the pseudo obs on phi instead on change in
### LoS directly

phiPseudo <- function(data, state.names, tra, cens.name, s = 0,
                       formula, timepoints, ncores = 1) {

    ## take care of the formula argument
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "id", "subset", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    Terms <- if (missing(data)) terms(formula)
    else terms(formula, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())

    ids <- unique(data$id)
    n <- length(ids)
    nt <- length(timepoints)

    ref <- as.matrix(predictPhi(clos(etm(data = data, state.names = state.names, tra = tra,
                                      cens.name = cens.name, s = 0, covariance = FALSE),
                                  aw = FALSE), timepoints)[, c("phi", "phi.case",
                                  "phi.control", "phiR")])

    ref <- apply(ref, 2, rep, n)
    psd <- matrix(0, nrow = n * nt, ncol = 6)
    
    temp <- parallel::mclapply(seq_along(ids), function(i) {
        as.matrix(predictPhi(clos(etm(data = data[!(data$id %in% ids[i]), ],
                                      state.names = state.names, tra = tra,
                                      cens.name = cens.name, s = 0, covariance = FALSE),
                                  aw = FALSE), timepoints)[, c("phi", "phi.case",
                                  "phi.control", "phiR")])
    }, mc.cores = ncores)
    ## } else {
    ##     temp <- lapply(seq_along(ids), function(i) {
    ##         as.matrix(predictPhi(clos(etm(data = data[!(data$id %in% ids[i]), ],
    ##                                    state.names = state.names, tra = tra,
    ##                                    cens.name = cens.name, s = 0, covariance = FALSE),
    ##                                aw = FALSE), timepoints)[, c("phi", "phi.case",
    ##                               "phi.control", "phiR")])
    ##     })
    ## }

    temp <- do.call(rbind, temp)

    for (i in seq_len(4)) {
        psd[, i + 2] <- n * ref[, i] - (n - 1) * temp[, i]
    }
    psd[, 1] <- as.vector(mapply(rep, ids, nt))
    psd[, 2] <- rep(timepoints, n)
    psd <- as.data.frame(psd)
    names(psd) <- c("id", "time", "ps.phi", "ps.phi.case",
                    "ps.phi.control", "ps.phiR")

    cov <- as.matrix(m[!duplicated(data$id), , drop = FALSE])
    cov <- matrix(mapply(rep, cov, nt), dim(psd)[1], dim(cov)[2])
    cov <- as.data.frame(cov)
    colnames(cov) <- attr(Terms, "term.labels")

    zzz <- list(pseudoData = data.frame(psd, cov),
                phi = data.frame(id = psd[, 1], ref, time = timepoints),
                ps = data.frame(id = psd[, 1], temp, time = timepoints),
                call = call)
    class(zzz) <- "phiPseudo"

    zzz
}

predictPhi <- function(object, timepoints) {
    if (!inherits(object, "clos.etm")) stop("gtfo")

    if (missing(timepoints)) stop("I want timepoints!!!")

    ## phi <- object$phi.case - object$phi.control
    
    ind <- findInterval(timepoints, object$time)
    tmp.case <- tmp.control <- numeric(length(timepoints))
    place <- which(ind != 0)
    tmp.case[place] <- object$phi.case[ind]
    tmp.control[place] <- object$phi.control[ind]

    data.frame(phi.case = tmp.case, phi.control = tmp.control,
               phi = tmp.case - tmp.control,
               phiR = tmp.case / tmp.control, time = timepoints)
}

    
