ci.land.vec <-
function (lambda, mu.hat, sig.sq.hat, nu, gamma.sq, ci.type = c("two-sided", 
    "lower", "upper"), conf.level = 0.95) 
{
    names.mu.hat <- names(mu.hat)
    ci.type <- match.arg(ci.type)
    if (any(abs(lambda) < .Machine$double.eps)) 
        stop("No values of 'lambda' can be 0")
    if (any(sig.sq.hat < .Machine$double.eps)) 
        stop("All values of 'sig.sq.hat' must be larger than 0")
    if (any(nu < 2) || any(nu != trunc(nu))) 
        stop("All values of 'nu' must be an integer greater than or equal to 2")
    if (any(gamma.sq < .Machine$double.eps)) 
        stop("All values of 'gamma.sq' must be larger than 0")
    if (length(conf.level) > 1 || conf.level < 0.5 || conf.level >= 
        1) 
        stop("'conf.level' must be a scalar between at least 50% and less than 100%")
    arg.mat <- cbind(lambda = as.vector(lambda), mu.hat = as.vector(mu.hat), 
        sig.sq.hat = as.vector(sig.sq.hat), nu = as.vector(nu), 
        gamma.sq = as.vector(gamma.sq))
    for (i in dimnames(arg.mat)[[2]]) assign(i, arg.mat[, i])
    n <- length(lambda)
    ret.obj <- matrix(0, n, 2)
    for (i in 1:n) ret.obj[i, ] <- ci.land(lambda[i], mu.hat[i], 
        sig.sq.hat[i], nu[i], gamma.sq[i], ci.type, conf.level)
    dimnames(ret.obj) <- list(names.mu.hat, c("LCL", "UCL"))
    attr(ret.obj, "ci.parameter") <- paste("mu +", lambda, "* (sigma^2)")
    attr(ret.obj, "ci.type") <- ci.type
    attr(ret.obj, "ci.method") <- "Land's (1971, 1975) Method"
    attr(ret.obj, "conf.level") <- conf.level
    ret.obj
}
