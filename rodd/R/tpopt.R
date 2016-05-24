tpopt <-  function(x, w = rep(1, length(x)) / length(x), eta, theta.fix, theta.var = NULL, p, x.lb = min(x), x.rb = max(x), opt = list())
{
    result <- topt(x, w, eta, theta.fix, theta.var, p, x.lb, x.rb, opt)
    class(result) <- "tpopt"
    result$call <- match.call()
    result
}

