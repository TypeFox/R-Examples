### To be used for competing endpoints
clos.cp <- function(x, tr.mat, aw, ratio) {
    dims <- dim(x$est)
    los <- matrix(rep(x$time, 3), ncol = 3, byrow = FALSE)
    phi2 <- matrix(data=c(x$time, rep(0, dims[3]), rep(0, dims[3])),
                   ncol=3, byrow=FALSE)
    phi3 <- matrix(data=c(x$time, rep(0, dims[3]), rep(0, dims[3])),
                   ncol=3, byrow=FALSE)
    ind.cens <- apply(x$n.event, 3, function(r) all(r == 0))
    tau <- max(x$time[ind.cens], x$time)
    
    out <- .C(los_cp,
              as.double(x$time),
              as.double(tr.mat),
              as.integer(dims[3]),
              as.integer(dims[1]),
              as.integer(dims[2]), 
              los1 = as.double(los[,2]),
              los0 = as.double(los[,3]),
              phi2case    = as.double(phi2[,2]),
              phi2control = as.double(phi2[,3]),
              phi3case    = as.double(phi3[,2]),
              phi3control = as.double(phi3[,3]),
              as.double(tau))
    
    los[, 2] <- out$los0
    los[, 3] <- out$los1
    phi2[, 3] <- out$phi2case; phi2[, 2] <- out$phi2control
    phi3[, 3] <- out$phi3case; phi3[, 2] <- out$phi3control
    indi <- apply(x$n.event, 3, function(x) {sum(x[1, ]) != 0})
    wait.times <- x$time[indi]
    wait.prob <- x$est["0", "0", ][indi]
    my.weights <- diff(c(0, 1 - wait.prob))
    
    pp <- x$n.risk[-1, ]
    ev.last <- apply(x$n.event[, , dims[3]], 1, sum)[1:2]
    pp <- rbind(pp, pp[nrow(pp), ] - ev.last)
    filtre <- pp[, 1] <= 0 | pp[, 2] <= 0
    
    tmp <- list(los, phi2, phi3)
    estimates <- lapply(tmp, function(z) {
        if (ratio) {
            ldiff <- z[, 3] / z[, 2]
        } else {
            ldiff <- z[, 3] - z[, 2]
        }
        ldiff[filtre] <- 0
        estimate <- matrix(ldiff[is.element(z[, 1], wait.times)], nrow = 1) %*%
            matrix(my.weights, ncol=1)
        estimate
    })
    
    e.phi.w1 <- e.phi.w23 <- my.weights1 <- my.weights23 <- NULL
    if (aw) {
        cif1 <- cumsum(c(1, x$est["0", "0", 1:(dims[3] - 1)]) * tr.mat[1, 2, ])
        my.weights1 <- diff(c(0, cif1[indi])) / cif1[length(cif1)]
        cif23 <- cumsum(c(1, x$est["0", "0", 1:(dims[3] - 1)]) *
                        (tr.mat[1, 3, ] + tr.mat[1, 4, ]))
        my.weights23 <- diff(c(0, cif23[indi])) / cif23[length(cif23)]
        weights.aw <- list(my.weights1, my.weights23)
        estimates.aw <- lapply(weights.aw, function(z) {
            ldiff <- los[, 3] - los[, 2]
            ldiff[filtre] <- 0
            estimate <- matrix(ldiff[is.element(los[, 1], wait.times)], nrow = 1) %*%
                matrix(z, ncol = 1)
            estimate
        })
        e.phi.w1 <- estimates.aw[[1]]
        e.phi.w23 <- estimates.aw[[2]]
    }
    
    res <- list(e.phi = estimates[[1]], phi.case = los[, 3],
                phi.control = los[, 2], e.phi2 = estimates[[2]],
                phi2.case = phi2[, 3], phi2.control = phi2[, 2],
                e.phi3 = estimates[[3]], phi3.case = phi3[, 3],
                phi3.control = phi3[, 2], weights = my.weights,
                w.time = wait.times, time = x$time, e.phi.weights.1 = e.phi.w1,
                e.phi.weights.other = e.phi.w23, weights.1 = my.weights1,
                weights.other = my.weights23)
    res
}


### To be used for single endpoint
clos.nocp <- function(x, tr.mat, aw, ratio) {
    dims <- dim(x$est)
    los <- matrix(rep(x$time, 3), ncol = 3, byrow = FALSE)
    tau <- max(x$time)
    
    out <- .C(los_nocp,
              as.double(x$time),
              as.double(tr.mat),
              as.integer(dims[3]),
              as.integer(dims[1]),
              as.integer(dims[2]), 
              los1 = as.double(los[,2]),
              los0 = as.double(los[,3]),
              as.double(tau))
    
    los[, 2] <- out$los0
    los[, 3] <- out$los1
    indi <- apply(x$n.event, 3, function(x) {sum(x[1, ]) != 0})
    wait.times <- x$time[indi]
    wait.prob <- x$est["0", "0", ][indi]

    pp <- x$n.risk[-1, ]
    ev.last <- apply(x$n.event[, , dims[3]], 1, sum)[1:2]
    pp <- rbind(pp, pp[nrow(pp), ] - ev.last)
    filtre <- pp[, 1] <= 0 | pp[, 2] <= 0

    if (ratio) {
        los.diff <- los[, 3] / los[, 2]
    } else {
        los.diff <- los[, 3] - los[, 2]
    }
    los.diff[filtre] <- 0
    my.weights <- diff(c(0, 1 - wait.prob))
    estimate <- matrix(los.diff[is.element(los[, 1], wait.times)], nrow = 1) %*%
        matrix(my.weights, ncol=1)
    
    e.phi.w1 <- e.phi.w2 <- my.weights1 <- my.weights2 <- NULL
    if (aw) {
        cif1 <- cumsum(c(1, x$est["0", "0", 1:(dims[3] - 1)]) * tr.mat[1, 2, ])
        my.weights1 <- diff(c(0, cif1[indi])) / cif1[length(cif1)]
        cif2 <- cumsum(c(1, x$est["0", "0", 1:(dims[3] - 1)]) * tr.mat[1, 3, ])
        my.weights2 <- diff(c(0, cif2[indi])) / cif2[length(cif2)]
        weights.aw <- list(my.weights1, my.weights2)
        estimates.aw <- lapply(weights.aw, function(z) {
            ldiff <- los[, 3] - los[, 2]
            ldiff[filtre] <- 0
            estimate <- matrix(ldiff[is.element(los[, 1], wait.times)], nrow = 1) %*%
                matrix(z, ncol = 1)
            estimate
        })
        e.phi.w1 <- estimates.aw[[1]]
        e.phi.w2 <- estimates.aw[[2]]
    }
    
    res <- list(e.phi = estimate[[1]], phi.case = los[, 3],
                phi.control = los[, 2], weights = my.weights,
                w.time = wait.times, time = x$time, e.phi.weights.1 = e.phi.w1,
                e.phi.weights.other = e.phi.w2, weights.1 = my.weights1,
                weights.other = my.weights2)
    res
}

    


clos <- function(x, aw = FALSE, ratio = FALSE) {
    if (!inherits(x, "etm")) {
        stop("'x' must be an 'etm' object")
    }
    if (is.null(x$delta.na)) {
        stop("Needs the increment of the Nelson-Aalen estimator")
    }
    absorb <- setdiff(levels(x$trans$to), levels(x$trans$from))
    transient <- unique(x$state.names[!(x$state.names %in% absorb)])
    if (!(length(transient) == 2 && length(absorb) %in% c(1, 2)))
        stop("The multistate model must have 2 transient states \n and 1 or 2 absorbing states")
    dims <- dim(x$est)
    comp.risk <- FALSE
    if (dims[1] == 4) comp.risk <- TRUE
    I <- diag(1, dims[1])
    tr.mat <- array(apply(x$delta.na, 3, "+", I), dim = dims)
    if (comp.risk) {
        res <- clos.cp(x, tr.mat, aw, ratio)
    }
    else res <- clos.nocp(x, tr.mat, aw, ratio)
    class(res) <- "clos.etm"
    res
}
