psi.lnorm.new <- function(x, ev.eta.f, sq.sigma.f, mu.f, eta.s, sq.var.s, theta.s)
{ 
    ev.eta.s <- eta.s(x, theta.s)
    sq.sigma.s <- log(1 + sq.var.s(x, theta.s) / (ev.eta.s * ev.eta.s))
    mu.s <- log(ev.eta.s) - 0.5 * sq.sigma.s

    temp <- mu.s - mu.f
                
    0.5 * (log(sq.sigma.f / sq.sigma.s) + (sq.sigma.s + temp * temp) / sq.sigma.f - 1)
}