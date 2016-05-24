KKT <- function(y, x, d, beta = NULL) {
    storage.mode(x) <- "double"
    ty <- as.double(y)
    tevent <- as.double(d)
    if (any(ty <= 0)) 
        stop("negative event times encountered;  not permitted for Cox family")
    nobs <- as.integer(length(ty))
    nvars <- as.integer(ncol(x))
    if (is.null(beta)) {
        beta <- double(0)
        nlam <- as.integer(1)
        nvars <- as.integer(0)
    } else {
        beta <- as.matrix(beta)
        nlam <- as.integer(ncol(beta))
    }
    fit <- .Fortran("KKT", nobs, nvars, as.double(x), ty, tevent, as.double(beta), 
        nlam, flog = double(nvars * nlam), PACKAGE = "fastcox")
    res <- matrix(fit$flog, nvars, nlam)
    res <- -res
}



KKTcheck <- function(y, x, d, alpha, p, lam, thr, beta) {
    fj <- KKT(y = y, x = x, d = d, beta = beta)
    B <- as.matrix(beta)
    ctr <- 0
    for (l in 1:length(lam)) {
        for (j in 1:p) {
            if (B[j, l] != 0) {
                AA <- fj[j, l] + lam[l] * (1 - alpha) * B[j, l] + alpha * lam[l] * 
                  sign(B[j, l])
                if (abs(AA) > thr) {
                  cat("violate b != 0", AA, "\n")
                  ctr <- ctr + 1
                  
                }
            } else {
                BB <- abs(fj[j, l]) - alpha * lam[l]
                if (BB > thr) {
                  cat("violate b = 0", BB, "\n")
                  ctr <- ctr + 1
                }
            }
        }
    }
    cat("# of violations", ctr, "\n")
    return(ctr)
}



OBJ <- function(y, x, d, beta = NULL) {
    storage.mode(x) <- "double"
    ty <- as.double(y)
    tevent <- as.double(d)
    if (any(ty <= 0)) 
        stop("negative event times encountered;  not permitted for Cox family")
    nobs <- as.integer(length(ty))
    nvars <- as.integer(ncol(x))
    if (is.null(beta)) {
        beta <- double(0)
        nlam <- as.integer(1)
        nvars <- as.integer(0)
    } else {
        beta <- as.matrix(beta)
        nlam <- as.integer(ncol(beta))
    }
    fit <- .Fortran("OBJ", nobs, nvars, as.double(x), ty, tevent, as.double(beta), 
        nlam, loss = double(nlam), PACKAGE = "fastcox")
    res <- fit$loss
    res
}



OBJcheck <- function(y, x, d, alpha, lam, beta) {
    o <- rep(0, length(lam))
    loss <- OBJ(y = y, x = x, d = d, beta = beta)
    B <- as.matrix(beta)
    ctr <- 0
    for (l in 1:length(lam)) {
        o[l] <- loss[l] + 0.5 * lam[l] * (1 - alpha) * crossprod(B[, l], B[, 
            l]) + alpha * lam[l] * sum(abs(B[, l]))
    }
    o
} 
