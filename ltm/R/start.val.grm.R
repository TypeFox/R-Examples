start.val.grm <-
function (start.val, data, weight, constrained, ncatg) {
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
            y <- data[, i]
            na.ind <- !is.na(y)
            y. <- y[na.ind]
            z. <- z[na.ind]
            weight. <- weight[na.ind]
            lev <- length(unique(y.))
            q <- lev - 1
            q1 <- lev %/% 2
            y1 <- (y. > q1)
            fit <- glm.fit(cbind(1, z.), y1, weight., family = binomial())
            coefs <- fit$coefficients
            spacing <- qlogis((1:q) / (q + 1))
            thets <- -coefs[1] + spacing - spacing[q1]
            out <- c(thets[1], log(diff(thets)), coefs[-1])
            names(out) <- NULL
            res[[i]] <- out
        }
        if (constrained)
            res[seq(1, p - 1)] <- lapply(res[seq(1, p - 1)], function (x) x[-length(x)])
        res
    } else {
        lapply(start.val, function (x) {
            nx <- length(x)
            c(x[1], log(diff(x[-nx])), x[nx])
        })
    }
}
