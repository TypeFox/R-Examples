KLopt.lnorm <- function(x, w = rep(1, length(x)) / length(x), eta, sq.var, theta.fix, theta.var = NULL, p, x.lb = min(x), x.rb = max(x), opt = list())
{
    result <- kl.opt.lnorm.main(x, w, eta, sq.var, theta.fix, theta.var, p, x.lb, x.rb, opt)
    class(result) <- "KLopt.lnorm"
    result$call <- match.call()
    result
}