tTestPower <-
function (n.or.n1, n2 = n.or.n1, delta.over.sigma = 0, alpha = 0.05, 
    sample.type = ifelse(!missing(n2), "two.sample", "one.sample"), 
    alternative = "two.sided", approx = FALSE) 
{
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    if (!is.vector(n.or.n1, mode = "numeric") || !is.vector(delta.over.sigma, 
        mode = "numeric") || !is.vector(alpha, mode = "numeric")) 
        stop("'n.or.n1', 'delta.over.sigma', and 'alpha' must be numeric vectors.")
    if (any(is.na(n.or.n1)) || any(is.na(delta.over.sigma))) 
        stop(paste("Missing (NA) and Undefined (Nan) values", 
            "are not allowed in 'n.or.n1' or 'delta.over.sigma'"))
    if (!all(is.finite(alpha))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'alpha'"))
    if (sample.type == "one.sample" && any(n.or.n1 < 2)) 
        stop("All values of 'n.or.n1' must be greater than or equal to 2.")
    if (any(alpha <= 0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be greater than 0 and less than 1.")
    if (sample.type == "two.sample") {
        if (!missing(n2)) {
            if (!is.vector(n2, mode = "numeric")) 
                stop("'n2' must be a numeric vector")
            if (any(is.na(n2))) 
                stop(paste("Missing (NA) and Undefined (Nan) values", 
                  "are not allowed in 'n2'"))
        }
        if (any(n.or.n1 < 1) || any(n2 < 1) || any((n.or.n1 + 
            n2) < 3)) 
            stop(paste("When sample.type='two.sample',", "all values of 'n.or.n1' and 'n2' must be", 
                "greater than or equal to 1 and they must sum", 
                "to at least 3"))
    }
    if (sample.type == "one.sample") {
        df <- n.or.n1 - 1
        ncp <- sqrt(n.or.n1) * delta.over.sigma
    }
    else {
        df <- n.or.n1 + n2 - 2
        ncp <- sqrt(1/(1/n.or.n1 + 1/n2)) * delta.over.sigma
    }
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
