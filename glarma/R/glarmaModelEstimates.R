glarmaModelEstimates <- function(object) {
    delta <- object$delta
    se <- sqrt(diag(object$cov))
    te <- delta/se
    pe <- 2 * (1 - pnorm(abs(te), 0, 1))
    estimates <- data.frame(delta, se, te, pe)
    rownames(estimates) <- names(delta)
    colnames(estimates) <- c("Estimate", "Std.Error", "z-ratio", "Pr(>|z|)")
    estimates
} 
