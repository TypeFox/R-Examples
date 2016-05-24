plot.KLopt.lnorm <- function(x, ...)
{
    res <- x
    curve(Psi.lnorm(x, res$p, res$eta, res$sq.var, res$theta.fix, res$theta.var), res$x.lb, res$x.rb, ylab = "Psi", n = 1001, type = "n")
    grid()
    curve(Psi.lnorm(x, res$p, res$eta, res$sq.var, res$theta.fix, res$theta.var), res$x.lb, res$x.rb, n = 1001, add = TRUE)
    points(res$x, Psi.lnorm(res$x, res$p, res$eta, res$sq.var, res$theta.fix, res$theta.var))
    abline(h = Psi.lnorm(res$x[1], res$p, res$eta, res$sq.var, res$theta.fix, res$theta.var), lty = 2)
}