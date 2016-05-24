start.val.tpm <-
function (start.val, data, type, constraint) {
    n <- nrow(data)
    p <- ncol(data)
    cmptStrVal <- is.null(start.val) || (start.val == "random" || (is.matrix(start.val) && length(start.val) != 3*p))
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
        coefs <- cbind(qlogis(seq(0.05, 0.15, length = p))[order(order(coefs[, 1], decreasing = TRUE))], coefs)
        coefs <- if (type == "rasch") c(coefs[, 1:2], abs(mean(coefs[, 3]))) else as.vector(coefs)
    } else {
        coefs <- start.val
        coefs[, 1] <- qlogis(coefs[, 1])
        coefs <- if (type == "rasch") c(coefs[, 1:2], abs(mean(coefs[, 3]))) else as.vector(coefs)
    }
    if (!is.null(constraint)) {
        if (type == "rasch" && any(ind <- constraint[, 2] == 3))
            coefs[-c((constraint[!ind, 2] - 1) * p + constraint[!ind, 1], length(coefs))]
        else
            coefs[-((constraint[, 2] - 1) * p + constraint[, 1])]
    } else
        coefs
}
