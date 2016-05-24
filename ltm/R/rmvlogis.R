rmvlogis <-
function (n, thetas, IRT = TRUE, link = c("logit", "probit"), 
        distr = c("normal", "logistic", "log-normal", "uniform"), z.vals = NULL) {
    if (!is.matrix(thetas) || !is.numeric(thetas))
        stop("'thetas' must be a numeric matrix with rows representing the items.\n")
    link <- match.arg(link)
    distr <- match.arg(distr)
    z <- if (is.null(z.vals) || length(z.vals) != n) {
        switch(distr, 
            "normal" = rnorm(n), 
            "logistic" = sqrt(3) / pi * rlogis(n),
            "log-normal" = (rlnorm(n) - exp(0.5)) / sqrt(exp(2) - exp(1)),
            "uniform" = runif(n, -3.5, 3.5) / sqrt(7^2/ 12))
    } else {
        c(z.vals)
    }
    p <- nrow(thetas)
    if (ncol(thetas) < 2 || ncol(thetas) > 3)
        stop("'thetas' must be either a two- or a three-column matrix.\n")
    betas <- thetas[, 1:2]
    eta <- if (IRT) {
        outer(z, betas[, 1], "-") * rep(betas[, 2], each = n)
    } else {
        cbind(1, z) %*% t(betas)
    }
    pr <- if (ncol(thetas) == 3) {
        cs <- thetas[, 3]
        if (any(cs < 0 | cs > 1))
            stop("some guessing parameters are either smaller than zero or greater than one.\n")
        cs.mat <- matrix(cs, n, p, TRUE)
        cs.mat + (1 - cs.mat) * if (link == "logit") plogis(eta) else pnorm(eta)
    } else {
        if (link == "logit") plogis(eta) else pnorm(eta)
    }
    X <- matrix(0, n, p)
    for (i in 1:p)
        X[, i] <- rbinom(n, 1, pr[, i])
    X
}
