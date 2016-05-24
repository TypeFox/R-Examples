"start.values" <-
function(a, b, X){
    old <- options(warn = (-1))
    on.exit(options(old))
    N <- round(ifelse(a == 0 | b == 1, 0.5 / (b - a), 1 / (b - a)))
    y <- N * (a + b) / 2
    y <- ifelse(a == 0, 0, y)
    y <- ifelse(b == 1, N, y)
    test <- try(glm(cbind(y, N - y) ~ X[, -1], family = binomial)$coef, silent = TRUE)
    start.par <- if(class(test) == "try-error") c(rnorm(ncol(X)), 1.) else c(test, 1.)
    names(start.par) <- NULL
    start.par
}

