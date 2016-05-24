mu <- function(x, eta.f, sq.var.f, theta.f)
{
    ev.eta.f <- eta.f(x, theta.f)
    sq.sigma.f <- log(1 + sq.var.f(x, theta.f) / (ev.eta.f * ev.eta.f))
    log(ev.eta.f) - 0.5 * sq.sigma.f
}