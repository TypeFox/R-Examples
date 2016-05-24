rmvordlogis <-
function (n, thetas, IRT = TRUE, model = c("gpcm", "grm"), link = c("logit", "probit"),
    distr = c("normal", "logistic", "log-normal", "uniform"), z.vals = NULL) {
    if (n < 0 || !is.list(thetas))
        stop("invalid arguments.\n")
    model <- match.arg(model)
    link <- match.arg(link)
    distr <- match.arg(distr)
    z <- if (is.null(z.vals) || length(z.vals) != n) {
        switch(distr, 
            "normal" = rnorm(n), 
            "logistic" = sqrt(3) / pi * rlogis(n),
            "log-normal" = (rlnorm(n) - exp(0.5)) / sqrt(exp(2) - exp(1)),
            "uniform" = runif(n, -3.5, 3.5) / sqrt(7^2/ 12))
    } else {
        z.vals
    }
    p <- length(thetas)
    ncatg <- sapply(thetas, length)
    probs <- if (model == "grm") {
        gammas <- lapply(thetas, function (x) {
            nx <- length(x)
            if (IRT)
                 cbind(plogis(x[nx] * (z - matrix(x[-nx], n, nx - 1, TRUE))), 1)
            else 
                cbind(plogis(matrix(x[-nx], n, nx - 1, TRUE) - x[nx] * z), 1)
        })
        lapply(gammas, function (x) {
            nc <- ncol(x)
            cbind(x[, 1], x[, 2:nc] - x[, 1:(nc - 1)])
        })
    } else {
        lapply(crf.GPCM(thetas, z, IRT), t)
    }
    X <- matrix(0, n, p)
    for (j in 1:p) {
        for (i in 1:n)
            X[i, j] <- sample(ncatg[j], 1, prob = probs[[j]][i, ])
    }
    X
}
