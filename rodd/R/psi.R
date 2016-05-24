psi <- function(x, eta.f, eta.s, theta.f, theta.s)
{
    temp <- (eta.f(x, theta.f) - eta.s(x, theta.s))
    temp * temp
}
