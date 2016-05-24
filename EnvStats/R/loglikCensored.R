loglikCensored <-
function (theta, x, censored, censoring.side = c("left", "right"), 
    distribution = "norm") 
{
    censoring.side <- match.arg(censoring.side)
    ret.val <- switch(censoring.side, right = sum(log(1 - do.call(paste("p", 
        distribution, sep = ""), args = c(list(q = x[censored]), 
        as.list(theta))))) + sum(log(do.call(paste("d", distribution, 
        sep = ""), args = c(list(x = x[!censored]), as.list(theta))))), 
        left = sum(log(do.call(paste("p", distribution, sep = ""), 
            args = c(list(q = x[censored]), as.list(theta))))) + 
            sum(log(do.call(paste("d", distribution, sep = ""), 
                args = c(list(x = x[!censored]), as.list(theta))))))
    ret.val
}
