epoisSinglyCensored.half.cen.level <-
function (x, censored, censoring.side, ci, ci.method = "normal.approx", 
    ci.type, conf.level, ci.sample.size = sum(!censored), pivot.statistic = c("z", 
        "t")) 
{
    x[censored] <- x[censored]/2
    lambda.hat <- c(lambda = mean(x))
    if (ci) {
        ci.method <- match.arg(ci.method)
        pivot.statistic <- match.arg(pivot.statistic)
        ci.obj <- ci.normal.approx(theta.hat = lambda.hat, sd.theta.hat = sqrt(lambda.hat/ci.sample.size), 
            n = ci.sample.size, df = ci.sample.size - 1, ci.type = ci.type, 
            alpha = 1 - conf.level, lb = 0, test.statistic = pivot.statistic)
        ci.obj$parameter <- "lambda"
        return(list(parameters = lambda.hat, ci.obj = ci.obj))
    }
    else return(list(parameters = lambda.hat))
}
