linearTrendTestPower <-
function (n, x = lapply(n, seq), slope.over.sigma = 0, alpha = 0.05, 
    alternative = "two.sided", approx = FALSE) 
{
    if (missing(n) && missing(x)) 
        stop("You must supply either 'n' or 'x'")
    if (!missing(x)) {
        if (is.vector(x) && !is.list(x)) 
            x <- list(x)
        if (!is.list(x) || !all(sapply(x, function(i) {
            is.vector(i, mode = "numeric")
        }))) 
            stop(paste("'x' must be either a numeric vector or", 
                "a list in which each component of 'x'", "is a numeric vector"))
        if (any(sapply(x, function(i) {
            !all(is.finite(i))
        }))) 
            stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                "Undefined (Nan) values are not allowed in", 
                "any components of 'x'"))
        n <- sapply(x, length)
        n.unique <- sapply(x, function(y) length(unique(y)))
        if (any(n < 3) || any(n.unique < 2)) 
            stop(paste("All components of 'x' must contain at least 3 elements", 
                "and at least 2 distinct values."))
    }
    else {
        if (!is.vector(n, mode = "numeric")) 
            stop("'n' must be a numeric vector")
        if (any(is.na(n))) 
            stop(paste("Missing (NA) and Undefined (Nan) values", 
                "are not allowed in 'n'"))
        if (any(n < 3)) 
            stop("All values of 'n' must be greater than or equal to 3.")
    }
    if (!is.vector(slope.over.sigma, mode = "numeric") || !is.vector(alpha, 
        mode = "numeric")) 
        stop("'slope.over.sigma', and 'alpha' must be numeric vectors.")
    if (any(is.na(slope.over.sigma))) 
        stop(paste("Missing (NA) and Undefined (Nan) values", 
            "are not allowed in 'slope.over.sigma'"))
    if (!all(is.finite(alpha))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'alpha'"))
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be between 0 and 1.")
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    df <- n - 2
    ncp <- sqrt(sapply(x, function(x) (length(x) - 1) * var(x))) * 
        slope.over.sigma
    if (approx) 
        power <- switch(alternative, less = pt(qt(alpha, df) - 
            ncp, df), greater = 1 - pt(qt(1 - alpha, df) - ncp, 
            df), two.sided = pt(qt(alpha/2, df) - ncp, df) + 
            1 - pt(qt(1 - alpha/2, df) - ncp, df))
    else power <- switch(alternative, less = pT(qt(alpha, df), 
        df = df, ncp = ncp), greater = 1 - pT(qt(1 - alpha, df), 
        df = df, ncp = ncp), two.sided = 1 - pf(qf(1 - alpha, 
        1, df), df1 = 1, df2 = df, ncp = ncp^2))
    power
}
