`cuspfit` <-
function (x, start = c(0, 1, 0, 1), ..., objf = cusp.logLike, 
    method = "L-BFGS-B", lower = if (method == "L-BFGS-B") c(-50, 
        -50, -50, 0) else -Inf, upper = c(50, 50, 50, Inf)) 
{
    s <- sd(x)
    m <- mean(x)
    x <- scale(x)
    fit <- optim(start, objf, x = x, ..., method = method, lower = lower, 
        upper = upper)
    fit$par[3] <- m + s * fit$par[3]
    fit$par[4] <- s * fit$par[4]
    fit
}

