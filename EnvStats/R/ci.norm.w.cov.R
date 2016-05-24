ci.norm.w.cov <-
function (muhat, sdhat, var.muhat, var.sdhat, cov.muhat.sdhat, 
    n, df, ci.type = c("two-sided", "lower", "upper"), alpha = 0.05) 
{
    ci.type <- match.arg(ci.type)
    ss <- sdhat^2
    gam.11 <- (var.muhat * n)/ss
    gam.22 <- (var.sdhat * n)/ss
    gam.12 <- (cov.muhat.sdhat * n)/ss
    if (ci.type == "two-sided") 
        u2 <- qnorm(alpha/2)^2
    else u2 <- qnorm(alpha)^2
    denom <- n - u2 * gam.22
    if (abs(denom) < .Machine$double.eps) 
        stop("Zero denominator")
    if (denom < .Machine$double.eps) 
        warning(paste("Due to the small sample size and large proportion censored,", 
            "this approximation for the CI yields reversed bounds."))
    con <- u2/denom
    A <- -gam.12 * con
    B <- gam.11 * con
    con <- sqrt(A^2 + B)
    if (ci.type == "two-sided") {
        z1 <- (-A) + con
        z2 <- A + con
        lcl <- muhat - z2 * sdhat
        ucl <- muhat + z1 * sdhat
    }
    else if (ci.type == "lower") {
        z2 <- A + con
        lcl <- muhat - z2 * sdhat
        ucl <- Inf
    }
    else {
        z1 <- (-A) + con
        lcl <- -Inf
        ucl <- muhat + z1 * sdhat
    }
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    ret.obj <- list(name = "Confidence", parameter = "mean", 
        limits = ci.limits, type = ci.type, method = "Normal Approximation w/ Covariance", 
        conf.level = 1 - alpha, sample.size = n, dof = df)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
