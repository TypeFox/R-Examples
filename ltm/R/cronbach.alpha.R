cronbach.alpha <-
function (data, standardized = FALSE, CI = FALSE, probs = c(0.025, 0.975), B = 1000, na.rm = FALSE) {
    if (!inherits(data, "matrix") && !inherits(data, "data.frame"))
        stop("'data' must be either a data.frame or a matrix.\n")
    n <- nrow(data)
    p <- ncol(data)
    nam <- deparse(substitute(data))
    if (p < 2)
        stop("'data' should have more than two columns.\n")
    data <- data.matrix(data)
    if (!na.rm && any(is.na(data))) {
        stop("missing values in 'data'.\n")
    } else {
        alpha <- if (!standardized) {
            VarTot <- var(rowSums(data[complete.cases(data), ]))
            VarInd <- sum(apply(data, 2, sd, na.rm = TRUE)^2)
            (p / (p - 1)) * (1 - (VarInd / VarTot))
        } else {
            mat <- cor(data, use = "complete.obs")
            ave.rho <- mean(mat[upper.tri(mat)])
            (p * ave.rho) / (1 + (p - 1) * ave.rho)
        }        
    }
    out <- list(alpha = alpha, n = n, p = p, standardized = standardized, name = nam)
    if (CI) {
        T.boot <- numeric(B)
        for (i in 1:B) {
            data.boot <- data[sample(1:n, replace = TRUE), ]
            T.boot[i] <- if (!standardized) {
                VarTot <- var(rowSums(data.boot[complete.cases(data.boot), ]))
                VarInd <- sum(apply(data.boot, 2, sd, na.rm = na.rm)^2)
                (p / (p - 1)) * (1 - (VarInd / VarTot))
            } else {
                mat <- if (na.rm) cor(data.boot, use = "complete.obs") else cor(data.boot)
                ave.rho <- mean(mat[upper.tri(mat)])
                (p * ave.rho) / (1 + (p - 1) * ave.rho)
            }
        }
        out$ci <- quantile(T.boot, probs = probs)
        out$probs <- probs
        out$B <- B
    }
    class(out) <- "cronbachAlpha"
    out
}
