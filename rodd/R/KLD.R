KLD <- function(x, w, eta.f, eta.s, sq.var.f, sq.var.s, theta.f, theta.s)
{   
    sum(w * psi.lnorm(x, eta.f, eta.s, sq.var.f, sq.var.s, theta.f, theta.s))
}