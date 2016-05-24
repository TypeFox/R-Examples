epoisMultiplyCensored.mle <-
function (x, censored, censoring.side, ci, ci.method = "profile.likelihood", 
    ci.type, conf.level, ci.sample.size = N - n.cen, pivot.statistic = "z") 
{
    N <- length(x)
    x.cen <- x[censored]
    n.cen <- length(x.cen)
    cen.levels <- sort(unique(x.cen))
    c.vec <- table(x.cen)
    x.bar <- mean(x[!censored])
    if (censoring.side == "left") 
        fcn <- function(lambda, x.bar, N, cen.levels, c.vec, 
            n.cen) {
            c1 <- c.vec/(N - n.cen)
            c2 <- dpois(cen.levels - 1, lambda)/ppois(cen.levels - 
                1, lambda)
            (x.bar - lambda * (1 + sum(c1 * c2)))^2
        }
    else fcn <- function(lambda, x.bar, N, cen.levels, c.vec, 
        n.cen) {
        c1 <- c.vec/(N - n.cen)
        c2 <- dpois(cen.levels, lambda)/(1 - ppois(cen.levels, 
            lambda))
        (x.bar - lambda * (1 - sum(c1 * c2)))^2
    }
    lambda.hat <- nlminb(start = x.bar, objective = fcn, lower = .Machine$double.eps, 
        x.bar = x.bar, N = N, cen.levels = cen.levels, c.vec = c.vec, 
        n.cen = n.cen)$par
    names(lambda.hat) <- "lambda"
    ret.list <- list(parameters = lambda.hat)
    if (ci) {
        ci.method <- match.arg(ci.method, c("normal.approx", 
            "profile.likelihood"))
        pivot.statistic <- match.arg(pivot.statistic, c("z", 
            "t"))
        n <- N - n.cen
        if (censoring.side == "left") {
            c1 <- ppois(cen.levels - 1, lambda.hat)
            c2 <- dpois(cen.levels - 1, lambda.hat)/c1
            d2.lnL.wrt.lambda <- (-n * x.bar)/lambda.hat^2 - 
                sum(c.vec * (dpois(cen.levels - 2, lambda.hat)/c1 - 
                  c2 + c2^2))
        }
        else {
            c1 <- 1 - ppois(cen.levels, lambda.hat)
            c2 <- dpois(cen.levels, lambda.hat)/c1
            d2.lnL.wrt.lambda <- (-n * x.bar)/lambda.hat^2 + 
                sum(c.vec * (dpois(cen.levels - 1, lambda.hat)/c1 - 
                  c2 - c2^2))
        }
        var.lambda.hat <- -1/d2.lnL.wrt.lambda
        var.cov.params <- var.lambda.hat
        names(var.cov.params) <- "lambda"
        ci.obj <- ci.normal.approx(theta.hat = lambda.hat, sd.theta.hat = sqrt(var.lambda.hat), 
            n = ci.sample.size, df = ci.sample.size - 1, ci.type = ci.type, 
            alpha = 1 - conf.level, lb = 0, test.statistic = pivot.statistic)
        ci.obj$parameter <- "lambda"
        if (ci.method == "profile.likelihood") {
            limits <- ci.obj$limits
            names(limits) <- NULL
            ci.obj <- ci.epoisCensored.profile.likelihood(x = x, 
                censored = censored, censoring.side = censoring.side, 
                lambda.mle = lambda.hat, ci.type = ci.type, conf.level = conf.level, 
                LCL.start = limits[1], UCL.start = limits[2])
        }
        ret.list <- c(ret.list, list(var.cov.params = var.cov.params, 
            ci.obj = ci.obj))
    }
    ret.list
}
