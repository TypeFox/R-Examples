KLD.new <- function(x, w, ev.eta.f, sq.sigma.f, mu.f, eta.s, sq.var.s, theta.s)
{   
    sum(w * psi.lnorm.new(x, ev.eta.f, sq.sigma.f, mu.f, eta.s, sq.var.s, theta.s))
}