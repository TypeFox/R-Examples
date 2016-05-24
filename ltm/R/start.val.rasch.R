start.val.rasch <-
function (start.val, data) {
    n <- nrow(data)
    p <- ncol(data)
    cmptStrVal <- is.null(start.val) || (start.val == "random" || (all(is.numeric(start.val)) && length(start.val) != p+1))
    randStrVal <- length(start.val) == 1 && start.val == "random"
    if (cmptStrVal) {
        rs <- as.vector(rowSums(data, na.rm = TRUE))
        len.uni <- length(unique(rs))
        rs <- factor(rs, labels = 1:len.uni)
        rs <- as.numeric(levels(rs))[as.integer(rs)]
        z <- cbind(1, seq(-3, 3, len = len.uni)[rs])
        if (randStrVal)
            z[, 2] <- rnorm(n)
        old <- options(warn = (2))
        on.exit(options(old))
        coefs <- matrix(0, p, 2)
        for (i in 1:p) {
            y <- data[, i]
            na.ind <- !is.na(y)
            y. <- y[na.ind]
            z. <- z[na.ind, ]
            fm <- try(glm.fit(z., y., family = binomial()), silent = TRUE)
            coefs[i, ] <- if (!inherits(fm, "try-error")) {
                fm$coef
            } else {
                c(0, 1)
            }
        }
        c(coefs[, 1], mean(coefs[, 2]))
    } else 
        start.val
}
