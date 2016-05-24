sans.cov <- function(i, object, trs.sep) {
    P <- object$est[trs.sep[i, 1], trs.sep[i, 2], ]
    time <- object$time
    n.event <- object$n.event[trs.sep[i, 1], trs.sep[i, 2], ]
    n.risk <- object$n.risk[, trs.sep[i, 1]]
    data.frame(P, time, n.risk, n.event)
}

avec.cov <- function(i, object, transfo, trs.sep, trs, level) {
    P <- object$est[trs.sep[i, 1], trs.sep[i, 2], ]
    time <- object$time
    n.event <- object$n.event[trs.sep[i, 1], trs.sep[i, 2], ]
    n.risk <- object$n.risk[, trs.sep[i, 1]]
    var <- object$cov[trs[[i]], trs[[i]], ]
    alpha <- qnorm(level + (1 - level) / 2)
    switch(transfo[i],
           "linear" = {
               lower <- P - alpha * sqrt(var)
               upper <- P + alpha * sqrt(var)
           },
           "log" = {
               lower <- exp(log(P) - alpha * sqrt(var) / P)
               upper <- exp(log(P) + alpha * sqrt(var) / P)
           },
           "cloglog" = {
               lower <- 1 - (1 - P)^(exp(alpha * (sqrt(var) / ((1 - P) * log(1 - P)))))
               upper <- 1 - (1 - P)^(exp(-alpha * (sqrt(var) / ((1 - P) * log(1 - P)))))
           },
           "log-log" = {
               lower <- P^(exp(-alpha * (sqrt(var) / (P * log(P)))))
               upper <- P^(exp(alpha * (sqrt(var) / (P * log(P)))))
           })
    lower <- pmax(lower, 0)
    upper <- pmin(upper, 1)
    data.frame(P, time, var, lower, upper, n.risk, n.event)
}
    

ci.transfo <- function(object, tr.choice, level = 0.95, transfo = "linear") {
    if (!inherits(object, "etm")) {
        stop ("'x' must be of class 'etm'")
    }
    lt <- length(tr.choice)
    trs <- tr.choice
    trs.sep <- lapply(trs, strsplit, split = " ")
    ## Fixing separation of states with names including a space
    for (i in seq_along(trs.sep)) {
        if (length(trs.sep[[i]][[1]]) == 2) {
            next
        } else {
            tt <- charmatch(trs.sep[[i]][[1]], object$state.names, nomatch = 0)
            trs.sep[[i]][[1]] <- object$state.names[tt]
        }
    }
    trs.sep <- matrix(unlist(trs.sep), length(trs.sep), 2, byrow = TRUE)
    if (length(transfo) != lt)
        transfo <- rep(transfo[1], lt)
    if (is.null(object$cov)) {
        res <- lapply(seq_len(lt), sans.cov, object = object, trs.sep = trs.sep)
    }
    else {
        res <- lapply(seq_len(lt), avec.cov, object = object, transfo = transfo,
                      trs.sep = trs.sep, trs = trs, level = level)
    }
    names(res) <- tr.choice
    res
}
        
