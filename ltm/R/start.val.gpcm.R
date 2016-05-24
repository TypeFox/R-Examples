start.val.gpcm <-
function (start.val, data, weight, constraint, ncatg, IRT.param) {
    n <- nrow(data)
    p <- ncol(data)
    computeStartVals <- function (start.val) {
        ind <- if (!is.null(start.val)) {
            if (!is.list(start.val) && start.val == "random")
                return(list(compute = TRUE, random = TRUE))
            if (!is.list(start.val) && length(start.val) != p) {
                warning("'start.val' not of proper type; random starting values are used instead.\n")
                TRUE
            } else if (!all(ncatg == sapply(start.val, length))) {
                warning("number of parameter in 'start.val' differ from the number of levels in 'data'; random starting values are used instead.\n")
                TRUE
            } else
                FALSE      
        } else 
            TRUE
        list(compute = ind, random = FALSE)
    }
    comp <- computeStartVals(start.val)
    if (comp$compute) {
        res <- vector("list", p)
        z <- if (comp$random) rnorm(n) else seq(-3, 3, length = n)[order(rowSums(data, na.rm = TRUE))]
        for (i in 1:p) {
            out <- try({
                y <- data[, i]
                na.ind <- !is.na(y)
                y. <- y[na.ind]
                z. <- z[na.ind]
                weight. <- weight[na.ind]
                lev <- length(unique(y.))
                q <- lev - 1
                q1 <- lev %/% 2
                y1 <- (y. > q1)
                fit <- glm.fit(cbind(1, z.), y1, family = binomial()) #glm.fit(cbind(1, z.), y1, weight., family = binomial())
                coefs <- fit$coefficients
                names(coefs) <- NULL
                spacing <- qlogis((1:q) / (q + 1))
                thets <- - coefs[1] + spacing - spacing[q1]
                if (IRT.param) c(thets, coefs[-1]) else c(rev(thets), coefs[-1])
            }, TRUE)
            res[[i]] <- if (!inherits(out, "try-error")) out else {
                lapply(1:p, function (i) {
                    ss <- seq(-0.1, 0.1, len = ncatg[i] - 1)
                    if (IRT.param) c(ss, 1) else c(rev(ss), 1)
                })
            }
        }
        if (constraint == "1PL") {
            res[seq(1, p - 1)] <- lapply(res[seq(1, p - 1)], function (x) x[-length(x)])
        }
        if (constraint == "rasch")
            res <- lapply(res, function (x) x[-length(x)])
        unlist(res)
    } else {
        if (constraint == "gpcm") {
            unlist(start.val)
        } else if (constraint == "1PL") {
            start.val[seq(1, p - 1)] <- lapply(start.val[seq(1, p - 1)], function (x) x[-length(x)])
            unlist(start.val)
        } else {
            start.val <- lapply(start.val, function (x) x[-length(x)])
            unlist(start.val)
        }
    }
}
