psi.lnorm <- function(x, eta.f, eta.s, sq.var.f, sq.var.s, theta.f, theta.s)
{
    ev.eta.f <- eta.f(x, theta.f) 
    ev.eta.s <- eta.s(x, theta.s)
    sq.sigma.f <- log(1 + sq.var.f(x, theta.f) / (ev.eta.f * ev.eta.f))
    sq.sigma.s <- log(1 + sq.var.s(x, theta.s) / (ev.eta.s * ev.eta.s))
    mu.f <- log(ev.eta.f) - 0.5 * sq.sigma.f
    mu.s <- log(ev.eta.s) - 0.5 * sq.sigma.s

    temp <- mu.s - mu.f
            
    0.5 * (log(sq.sigma.f / sq.sigma.s) + (sq.sigma.s + temp * temp) / sq.sigma.f - 1)
}