information <-
function (object, range, items = NULL, ...) {
    if (!class(object) %in% c("grm", "gpcm", "ltm", "rasch", "tpm"))
        stop("'object' must inherit from either class 'grm', class 'gpcm', class 'ltm', class 'rasch' or class 'tpm'.\n")
    p <- ncol(object$X)
    itms <- if (!is.null(items)) {
        if (!is.numeric(items) && length(items) > p)
            stop("'items' should be a numeric vector of maximum length ", p, ".\n")
        if (any(!items %in% 1:p))
            stop("'items' should contain numbers in: ", paste(1:p, collapse = ", "), " indicating the items.\n")
        items
    } else 
        1:p
    if (class(object) == "ltm" && (object$ltst$factors > 1 | any(unlist(object$ltst[c("inter", "quad.z1", "quad.z2")]))))
        stop("Information is currently computed only for the two-parameter logistic model.\n")
    f <- function (z) {
        switch(class(object),
            "grm" = rowSums(infoprobs(object$coefficients, z)[, itms, drop = FALSE]),
            "gpcm" = rowSums(infoGPCM(object$coefficients, z, object$IRT.param)[, itms, drop = FALSE]),
            "ltm" = { betas <- object$coefficients; Z <- cbind(1, z)
                mat <- t(t(plogis(Z %*% t(betas)) * (1 - plogis(Z %*% t(betas)))) * betas[, 2]^2)
                rowSums(mat[, itms, drop = FALSE])
                },
            "rasch" = { betas <- object$coefficients; Z <- cbind(1, z)
                mat <- betas[1, 2]^2 * plogis(Z %*% t(betas)) * (1 - plogis(Z %*% t(betas)))
                rowSums(mat[, itms, drop = FALSE])
                },
            "tpm" = { thetas <- object$coefficients; Z <- cbind(1, z)
                betas <- thetas[, 2:3]
                cs <- plogis(thetas[, 1]) * object$max.guessing
                pi. <- plogis(Z %*% t(betas))
                cs <- matrix(cs, nrow(Z), p, TRUE)
                pi <- cs + (1 - cs) * pi.
                pqr <- pi * (1 - pi) * (pi. / pi)^2
                mat <- t(t(pqr) * betas[, 2]^2)
                rowSums(mat[, itms, drop = FALSE])
                })
    }
    I0 <- integrate(f, -10, 10, ...)$value
    I1 <- integrate(f, range[1], range[2], ...)$value
    out <- list("InfoRange" = I1, "InfoTotal" = I0, "PropRange" = I1 / I0, range = range, items = items, 
                call = object$call)
    class(out) <- "information"
    out
}
