dhsvm <- function(v, delta) {
    r <- v[1]
    if (r > 1) 
        dl <- 0 else if (r <= (1 - delta)) 
        dl <- -1 else dl <- (r - 1)/delta
    dl
}


dlogit <- function(r, delta) {
    dl <- -1/(1 + exp(r))
    dl
}

dsqsvm <- function(r, delta) {
    dl <- -2 * ifelse((1 - r) > 0, (1 - r), 0)
    dl
}

dls <- function(r, delta) {
    dl <- -r
}

margin <- function(b0, beta, y, x, delta, loss = c("ls", "logit", 
    "sqsvm", "hsvm")) {
    loss <- match.arg(loss)
    nobs <- nrow(x)
    b0MAT <- matrix(rep(b0, nobs), nrow = nobs, byrow = TRUE)
    link <- x %*% beta + b0MAT
    if (loss %in% c("logit", "sqsvm", "hsvm")) {
        r <- y * link
    } else r <- y - link
    fun <- paste("d", loss, sep = "")
    dMat <- apply(r, c(1, 2), eval(fun), delta = delta)
    if (loss %in% c("logit", "sqsvm", "hsvm")) {
        yxdMat <- t(x) %*% (dMat * y)/nobs
    } else yxdMat <- t(x) %*% dMat/nobs
    yxdMat
}


KKT <- function(b0, beta, y, x, lambda, pf, group, thr, delta, loss = c("ls", 
    "logit", "sqsvm", "hsvm")) {
    loss <- match.arg(loss)
    bn <- as.integer(max(group))
    dl <- margin(b0, beta, y, x, delta, loss)
    B <- matrix(NA, ncol = length(lambda))
    ctr <- 0
    for (l in 1:length(lambda)) {
        for (g in 1:bn) {
            ind <- (group == g)
            dl_norm <- sqrt(crossprod(dl[ind, l], dl[ind, l]))
            b_norm <- sqrt(crossprod(beta[ind, l], beta[ind, l]))
            if (b_norm != 0) {
                AA <- dl[ind, l] + beta[ind, l] * lambda[l] * pf[g]/b_norm
                if (abs(sum(AA)) >= thr) {
                  cat("violate at b != 0", abs(sum(AA)), "\n")
                  ctr <- ctr + 1
                }
            } else {
                BB <- dl_norm - pf[g] * lambda[l]
                if (BB > thr) {
                  cat("violate at b = 0", BB, "\n")
                  ctr <- ctr + 1
                }
            }
        }
    }
    cat("# of violations", ctr/length(lambda), "\n")
    return(ctr/length(lambda))
} 
