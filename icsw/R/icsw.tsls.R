icsw.tsls <- function(
    D, X, Y, Z, W,
    weights = NULL, R = 0,
    estimand = c("ATE", "ATT"),
    min.prob.quantile = NULL,
    min.prob = NULL, ...
) {
    n.obs <- length(Y)
    if (is.null(weights))
        weights <- rep(1, n.obs)
    
    fitted.model <- icsw.tsls.fit(
        D, X, Y, Z, W,
        weights = weights,
        estimand = estimand,
        min.prob = min.prob,
        min.prob.quantile = min.prob.quantile,
        ...
    )
    
    # bootstrap replicates
    if (R > 0) {
        coefs.boot <- matrix(
            nrow = 0,
            ncol = length(fitted.model$coefficients)
        )
        
        for (r in 1:R) {
            weights.r <- rexp(n.obs, 1) # Bayesian bootstrap
            coefs.r <- icsw.tsls.fit(
                D, X, Y, Z, W,
                weights = weights * weights.r,
                estimand = estimand,
                min.prob = min.prob,
                min.prob.quantile = min.prob.quantile,
                ...
            )$coefficients
            coefs.boot <- rbind(coefs.boot, coefs.r)
        }
        rownames(coefs.boot) <- 1:R
        return(c(
            fitted.model,
            list(
                coefs.boot = coefs.boot,
                coefs.se.boot = apply(coefs.boot, 2, sd)
            )
        ))
    } else {
        return(fitted.model)
    }
}