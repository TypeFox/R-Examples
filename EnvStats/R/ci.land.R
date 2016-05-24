ci.land <-
function (lambda, mu.hat, sig.sq.hat, n, nu, gamma.sq, ci.type = c("two-sided", 
    "lower", "upper"), conf.level = 0.95) 
{
    ci.type <- match.arg(ci.type)
    if (any(length.list(lambda, mu.hat, sig.sq.hat, nu, gamma.sq, 
        ci.type, conf.level) > 1)) 
        stop(paste("'lambda', 'mu.hat', 'sig.sq.hat', 'nu', 'ci.type', and 'conf.level'", 
            "must be vectors of length 1"))
    if (abs(lambda) < .Machine$double.eps) 
        stop("'lambda' cannot be 0")
    if (sig.sq.hat < .Machine$double.eps) 
        stop("'sig.sq.hat' must be larger than 0")
    if (nu < 2 || nu != trunc(nu)) 
        stop("'nu' must be an integer greater than or equal to 2")
    if (gamma.sq < .Machine$double.eps) 
        stop("'gamma.sq' must be larger than 0")
    if (conf.level < 0.5 || conf.level >= 1) 
        stop("'conf.level' must be at least 50% and less than 100%")
    k <- (nu + 1)/(2 * lambda * gamma.sq)
    S <- sqrt((2 * lambda * sig.sq.hat)/k)
    if (ci.type == "two-sided") {
        alpha <- (1 - conf.level)/2
        lcl <- mu.hat + lambda * sig.sq.hat + ((k * S)/sqrt(nu)) * 
            lands.C(S, nu, alpha)
        ucl <- mu.hat + lambda * sig.sq.hat + ((k * S)/sqrt(nu)) * 
            lands.C(S, nu, 1 - alpha)
    }
    else {
        if (ci.type == "lower") {
            lcl <- mu.hat + lambda * sig.sq.hat + ((k * S)/sqrt(nu)) * 
                lands.C(S, nu, 1 - conf.level)
            ucl <- Inf
        }
        else {
            lcl <- -Inf
            ucl <- mu.hat + lambda * sig.sq.hat + ((k * S)/sqrt(nu)) * 
                lands.C(S, nu, conf.level)
        }
    }
    ci.limits <- c(lcl, ucl)
    names(ci.limits) <- c("LCL", "UCL")
    ret.obj <- list(name = "Confidence", parameter = paste("mu +", 
        lambda, "* (sigma^2)"), limits = ci.limits, type = ci.type, 
        method = "Land", conf.level = conf.level, sample.size = n, 
        dof = nu)
    oldClass(ret.obj) <- "intervalEstimate"
    ret.obj
}
