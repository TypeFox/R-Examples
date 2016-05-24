start.val.ltm <-
function (start.val, data, factors, formula) {
    n <- nrow(data)
    p <- ncol(data)
    cf <- paste(formula[3])
    #form <- paste("y ~ ", paste("z", 1:factors, collapse = " + ", sep = ""), " + ", cf)
    form <- paste("y ~ ", cf)
    form <- as.formula(form)
    q. <- length(attr(terms(form), "term.labels")) + 1
    cmptStrVal <- is.null(start.val) || (start.val == "random" || (all(is.numeric(start.val)) && length(start.val) != p*q.))
    randStrVal <- length(start.val) == 1 && start.val == "random"
    if (cmptStrVal) {
        if (randStrVal) {
            Z <- data.frame(z1 = rnorm(n))
            if (factors > 1)
                Z$z2 <- rnorm(n)
        } else {
            rs <- as.vector(rowSums(data, na.rm = TRUE))
            len.uni <- length(unique(rs))
            rs <- factor(rs, labels = 1:len.uni)
            rs <- as.numeric(levels(rs))[as.integer(rs)]
            Z <- data.frame(z1 = seq(-3, 3, len = len.uni)[rs])
            if (factors > 1)
                Z$z2 <- seq(3, -3, len = n)
        }
        old <- options(warn = (2))
        on.exit(options(old))
        coefs <- matrix(0, p, q.)
        for (i in 1:p) {
            Z$y <- data[, i]
            fm <- try(glm(form, family = binomial(), data = Z), silent = TRUE)
            coefs[i, ] <- if (!inherits(fm, "try-error")) {
                fm$coef
            } else {
                c(0, rep(1, q. - 1))
            }
        }
        dimnames(coefs) <- NULL
        coefs
    } else 
        start.val
}
